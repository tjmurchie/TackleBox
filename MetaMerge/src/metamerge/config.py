"""Configuration loading and deep-merging for MetaMerge.

User config files are YAML documents that override individual keys in the
package defaults.  Only keys that differ from the defaults need to be listed;
unspecified keys retain their default values.

Example minimal config override file::

    thresholds:
      damage_min: 0.02
      strong_count_min_libraries: 3

    io:
      output_prefix: my_project

See ``config/defaults.yaml`` for the full list of available keys and their
default values, with inline documentation.
"""

from __future__ import annotations

import copy
from pathlib import Path
from typing import Any, Mapping

import yaml

from .defaults import DEFAULT_CONFIG


def deep_update(base: dict, updates: Mapping[str, Any]) -> dict:
    """Recursively merge ``updates`` into ``base``, returning a new dict.

    Nested dicts are merged recursively so that a partial override (e.g.,
    overriding only ``thresholds.damage_min``) does not erase sibling keys.
    Non-dict values (lists, scalars) are replaced wholesale.

    Args:
        base: The base dictionary (e.g., DEFAULT_CONFIG).
        updates: A mapping of overrides to apply.

    Returns:
        A new dictionary with ``updates`` applied on top of ``base``.
    """
    merged = copy.deepcopy(base)
    for key, value in updates.items():
        if isinstance(value, Mapping) and isinstance(merged.get(key), Mapping):
            merged[key] = deep_update(merged[key], value)
        else:
            merged[key] = value
    return merged


def load_config(config_path: str | None) -> dict:
    """Load the effective configuration for a MetaMerge run.

    Starts from the package defaults and deep-merges any user-supplied YAML
    config on top.  Returns the merged config dict.

    Args:
        config_path: Path to a user YAML config file, or ``None`` to use only
            the package defaults.

    Returns:
        Complete configuration dict ready to pass to workflow functions.

    Raises:
        FileNotFoundError: If ``config_path`` is set but the file does not exist.
        yaml.YAMLError: If the config file contains invalid YAML.
    """
    config = copy.deepcopy(DEFAULT_CONFIG)
    if config_path:
        path = Path(config_path)
        with path.open("r", encoding="utf-8") as handle:
            user = yaml.safe_load(handle) or {}
        config = deep_update(config, user)
    return config
