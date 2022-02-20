from covigator.database.queries import Queries
import plotly.graph_objects as go


PLOTLY_CONFIG = {
    'displaylogo': False,
    'displayModeBar': False,
    'modeBarButtonsToRemove': ['zoom', 'pan', 'select', 'zoomIn', 'zoomOut', 'autoScale', 'resetScale', 'lasso2d']
}
MARGIN = go.layout.Margin(l=0, r=0, b=0, t=30)
TEMPLATE = "plotly_white"

STYLES_STRIPPED = [{
                'if': {'row_index': 'odd'},
                'backgroundColor': 'rgb(248, 248, 248)'
            }]
STYLE_HEADER = {
    'backgroundColor': 'white',
    'fontWeight': 'bold'
}
STYLE_CELL = {
    'padding': '5px',
    'maxWidth': '100px'
}


class Figures:

    def __init__(self, queries: Queries):
        self.queries = queries
