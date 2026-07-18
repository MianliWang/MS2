"""Deterministic compound-and-batch isolated split helpers.

Rows are never shuffled or randomly divided.  A validation/final batch is held
out in full, and compounds occurring in that held-out set are purged from the
corresponding training/development set.
"""

from __future__ import annotations

import math
from collections.abc import Iterable, Sequence
from dataclasses import dataclass

from .features import FEATURE_NAMES, feature_vector


@dataclass(frozen=True)
class GroupedExample:
    """One eligible binary-truth example with immutable grouping identifiers."""

    row_id: str
    compound_id: str
    batch_id: str
    label: int
    features: tuple[float, ...]

    def __post_init__(self) -> None:
        if not str(self.row_id).strip():
            raise ValueError("row_id must be non-empty.")
        if not str(self.compound_id).strip():
            raise ValueError("compound_id must be non-empty.")
        if not str(self.batch_id).strip():
            raise ValueError("batch_id must be non-empty.")
        if self.label not in {0, 1}:
            raise ValueError("GroupedExample label must be 0 or 1.")
        normalized = feature_vector(self.features)
        if len(normalized) != len(FEATURE_NAMES) or not all(map(math.isfinite, normalized)):
            raise ValueError("GroupedExample features must be the six finite fixed features.")
        object.__setattr__(self, "features", normalized)


@dataclass(frozen=True)
class FoldSplit:
    """A validation batch with compound-overlap rows purged from training."""

    name: str
    train_indices: tuple[int, ...]
    validation_indices: tuple[int, ...]
    purged_indices: tuple[int, ...]

    def as_dict(self) -> dict[str, str | list[int]]:
        return {
            "name": self.name,
            "train_indices": list(self.train_indices),
            "validation_indices": list(self.validation_indices),
            "purged_indices": list(self.purged_indices),
        }


@dataclass(frozen=True)
class LockedBatchSplit:
    """One independent final batch and the leakage-purged development pool."""

    final_batch: str
    development_indices: tuple[int, ...]
    final_indices: tuple[int, ...]
    purged_indices: tuple[int, ...]

    def as_dict(self) -> dict[str, str | list[int]]:
        return {
            "final_batch": self.final_batch,
            "development_indices": list(self.development_indices),
            "final_indices": list(self.final_indices),
            "purged_indices": list(self.purged_indices),
        }


def build_locked_batch_split(
    examples: Sequence[GroupedExample],
    final_batch: str,
) -> LockedBatchSplit:
    """Hold out ``final_batch`` and purge its compounds from development."""

    _validate_examples(examples)
    final_name = str(final_batch).strip()
    if not final_name:
        raise ValueError("final_batch must be non-empty.")
    final_indices = tuple(
        index for index, example in enumerate(examples) if example.batch_id == final_name
    )
    if not final_indices:
        raise ValueError(f"Locked final batch {final_name!r} is absent from eligible truth data.")
    final_index_set = set(final_indices)
    final_compounds = {examples[index].compound_id for index in final_indices}
    development: list[int] = []
    purged: list[int] = []
    for index, example in enumerate(examples):
        if index in final_index_set:
            continue
        if example.compound_id in final_compounds:
            purged.append(index)
        else:
            development.append(index)
    if not development:
        raise ValueError("No development examples remain after locked-batch compound purge.")
    assert_group_isolation(examples, development, final_indices)
    return LockedBatchSplit(
        final_batch=final_name,
        development_indices=tuple(development),
        final_indices=final_indices,
        purged_indices=tuple(purged),
    )


def build_leave_one_batch_out_folds(
    examples: Sequence[GroupedExample],
    indices: Iterable[int] | None = None,
    name_prefix: str = "outer",
) -> tuple[FoldSplit, ...]:
    """Build deterministic batch folds with compound isolation in every fold.

    Calling this function on an outer fold's ``train_indices`` creates the
    corresponding nested grouped-CV inner folds.  There is intentionally no
    random-row fallback when too few batches remain.
    """

    _validate_examples(examples)
    pool = _normalize_indices(len(examples), indices)
    batches = sorted({examples[index].batch_id for index in pool})
    if len(batches) < 2:
        raise ValueError("Grouped cross-validation requires at least two eligible batches.")
    folds: list[FoldSplit] = []
    for validation_batch in batches:
        validation = tuple(index for index in pool if examples[index].batch_id == validation_batch)
        validation_compounds = {examples[index].compound_id for index in validation}
        candidate_train = [index for index in pool if examples[index].batch_id != validation_batch]
        purged = tuple(
            index
            for index in candidate_train
            if examples[index].compound_id in validation_compounds
        )
        purged_set = set(purged)
        train = tuple(index for index in candidate_train if index not in purged_set)
        if not train:
            raise ValueError(
                f"Fold {validation_batch!r} has no training rows after compound purge."
            )
        assert_group_isolation(examples, train, validation)
        folds.append(
            FoldSplit(
                name=f"{name_prefix}:{validation_batch}",
                train_indices=train,
                validation_indices=validation,
                purged_indices=purged,
            )
        )
    return tuple(folds)


def assert_group_isolation(
    examples: Sequence[GroupedExample],
    train_indices: Iterable[int],
    validation_indices: Iterable[int],
) -> None:
    """Raise if either compound IDs or batch IDs cross a split boundary."""

    train = _normalize_indices(len(examples), train_indices)
    validation = _normalize_indices(len(examples), validation_indices)
    overlap = set(train) & set(validation)
    if overlap:
        raise ValueError(f"Split reuses row indices across roles: {sorted(overlap)}.")
    train_compounds = {examples[index].compound_id for index in train}
    validation_compounds = {examples[index].compound_id for index in validation}
    compound_overlap = train_compounds & validation_compounds
    if compound_overlap:
        raise ValueError("Compound leakage across split: " + ", ".join(sorted(compound_overlap)))
    train_batches = {examples[index].batch_id for index in train}
    validation_batches = {examples[index].batch_id for index in validation}
    batch_overlap = train_batches & validation_batches
    if batch_overlap:
        raise ValueError("Batch leakage across split: " + ", ".join(sorted(batch_overlap)))


def validate_binary_fold(
    examples: Sequence[GroupedExample],
    split: FoldSplit | LockedBatchSplit,
) -> None:
    """Require positive and negative truth on both sides of a prospective test."""

    if isinstance(split, FoldSplit):
        training = split.train_indices
        evaluation = split.validation_indices
        evaluation_name = split.name
    else:
        training = split.development_indices
        evaluation = split.final_indices
        evaluation_name = f"final:{split.final_batch}"
    for role, indices in (("training", training), ("evaluation", evaluation)):
        labels = {examples[index].label for index in indices}
        if labels != {0, 1}:
            raise ValueError(
                f"{evaluation_name} {role} requires both truth classes; observed {sorted(labels)}."
            )


def _validate_examples(examples: Sequence[GroupedExample]) -> None:
    if not examples:
        raise ValueError("At least one eligible grouped example is required.")
    row_ids = [example.row_id for example in examples]
    if len(row_ids) != len(set(row_ids)):
        raise ValueError("Grouped example row_id values must be unique.")


def _normalize_indices(length: int, indices: Iterable[int] | None) -> tuple[int, ...]:
    if indices is None:
        return tuple(range(length))
    normalized = tuple(sorted({int(index) for index in indices}))
    if not normalized:
        raise ValueError("A split role cannot be empty.")
    if normalized[0] < 0 or normalized[-1] >= length:
        raise IndexError("Split index is outside the example sequence.")
    return normalized


__all__ = [
    "FoldSplit",
    "GroupedExample",
    "LockedBatchSplit",
    "assert_group_isolation",
    "build_leave_one_batch_out_folds",
    "build_locked_batch_split",
    "validate_binary_fold",
]
