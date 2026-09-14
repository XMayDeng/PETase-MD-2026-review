"""Missing measurements are errors, never implicit zero-valued observations."""


def require_complete(table):
    if table.isna().any().any():
        raise ValueError('Missing matrix entries: source data must explicitly encode every valid zero')
    return table
