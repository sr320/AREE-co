def _escape_cell(value):
    return str(value).replace("|", "\\|").replace("\r\n", " ").replace("\n", " ").replace("\r", " ")


def rows_to_markdown(columns, rows):
    if not rows:
        return "_No rows._"
    lines = [
        "| " + " | ".join(_escape_cell(col) for col in columns) + " |",
        "| " + " | ".join(["---"] * len(columns)) + " |",
    ]
    for row in rows:
        lines.append("| " + " | ".join(_escape_cell(value) for value in row) + " |")
    return "\n".join(lines)


def dataframe_to_markdown(df):
    columns = list(df.columns)
    return rows_to_markdown(columns, [[row[column] for column in columns] for _, row in df.iterrows()])
