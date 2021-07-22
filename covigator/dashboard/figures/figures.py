from covigator.database.queries import Queries
import plotly.graph_objects as go


PLOTLY_CONFIG = {
    'displaylogo': False,
    'displayModeBar': False,
    'modeBarButtonsToRemove': ['zoom', 'pan', 'select', 'zoomIn', 'zoomOut', 'autoScale', 'resetScale', 'lasso2d']
}
MARGIN = go.layout.Margin(l=0, r=0, b=0, t=30)
TEMPLATE = "plotly_white"


class Figures:

    def __init__(self, queries: Queries):
        self.queries = queries
