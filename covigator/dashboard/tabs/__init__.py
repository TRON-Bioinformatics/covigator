from datetime import datetime

import dash_html_components as html

TAB_STYLE = {
    'borderBottom': '1px solid #d6d6d6',
    'padding': '6px',
    'fontWeight': 'bold'
}
TAB_SELECTED_STYLE = {
    'borderTop': '1px solid #d6d6d6',
    'borderBottom': '1px solid #d6d6d6',
    'backgroundColor': '#119DFF',
    'color': 'white',
    'padding': '6px'
}
MONTH_PATTERN = "%Y-%m"
MISSING_VALUE = "-"
COLOR_OVERVIEW_MINI_CONTAINER = "#82B1FF"


def print_date(date: datetime.date):
    return str(date) if date is not None else MISSING_VALUE


def print_number(value):
    return '{:,}'.format(value)


def get_mini_container(title, value, color="#f9f9f9"):
    return html.Div(
        children=[
            html.P(html.B(title)),
            html.P(value)
        ],
        className="mini_container",
        style={"background-color": color, "font-size": 16}
    )
