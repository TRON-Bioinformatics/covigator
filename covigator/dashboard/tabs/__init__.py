from datetime import datetime
import dash_bootstrap_components as dbc

MONTH_PATTERN = "%Y-%m"
MISSING_VALUE = "-"


def print_date(date: datetime.date):
    return str(date) if date is not None else MISSING_VALUE


def print_number(value):
    return '{:,}'.format(value)


def get_mini_container(title, value):
    return dbc.Button(
        "{title}: {value}".format(title=title, value=value),
        outline=True,
        color="dark",
        className="me-1",
        style={"margin-right": "30px", "margin-left": "30px", 'font-size': '110%'})
