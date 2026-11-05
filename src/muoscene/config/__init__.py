"""Configuration schema and loading.

The pydantic schema is the single source of truth for every option and doubles as the
configuration reference in the documentation. The loader reads YAML, resolves `include`,
makes paths relative to the including file, fills defaults and produces the fully resolved
configuration stored in each run folder. The registry maps type names (detector types,
flux models, file formats) to implementations.

This package validates only; objects are built by each module's `from_config`.

Layer 0: may import only `muoscene` modules of layer <= 0 (see `muoscene._layers`).
"""
