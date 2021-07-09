from flask import Flask
from covigator.configuration import Configuration
from covigator.dashboard.dashboard import Dashboard
from logzero import logger

dashboard = Dashboard(config=Configuration())
app = dashboard.get_application()
server = Flask(__name__)
app.init_app(server)

if __name__ == '__main__':
    logger.info("Starting...")
    server.run()
    logger.info("Running...")
