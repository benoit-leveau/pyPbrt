"""Logging object used in the entire application."""


import logging


FORMAT = '%(asctime)-15s  %(message)s'
logger = None

def init_logger():
    global logger
    
    # setting the format    
    logging.basicConfig(format=FORMAT)

    # constructs the logger
    logger = logging.getLogger('pypbrt')

# init
init_logger()
