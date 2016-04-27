
import leip
import logging

def dispatch():
    """
    Run the nfj app
    """
    app.run()

logging.getLogger('nfj').setLevel(logging.INFO)

app = leip.app(name='nfj')
app.discover(globals())
