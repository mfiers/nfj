
import leip

@leip.arg('name', help='Say hello to')
@leip.command
def hello_world(app, args):
    print("{} {}".format(app.conf['message'], args.name))
