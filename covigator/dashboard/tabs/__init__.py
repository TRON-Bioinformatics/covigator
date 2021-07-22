from datetime import datetime
import dash_bootstrap_components as dbc

MONTH_PATTERN = "%Y-%m"
MISSING_VALUE = "-"


def print_date(date: datetime.date):
    return str(date) if date is not None else MISSING_VALUE


def print_number(value):
    return '{:,}'.format(value)


def get_mini_container(title, value, font_size=20):
    return dbc.Button(
        [title + " ", dbc.Badge(value, color="light", className="ml-1")],
        color="warning", size='lg',  style={"margin-left": "20px", "font-size": font_size})
