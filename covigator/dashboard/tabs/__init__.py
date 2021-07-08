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


def get_mini_container(title, value):
    return html.Div(
        children=[
            html.H6(title),
            html.H6(value)
        ],
        className="mini_container",
    )
