import sys, getopt
import matplotlib.pyplot as plt
import getSpecies as gs

def plot(tags, tmin, tmax, Tmin, Tmax, scale):
    """ Plotting routine """
    fields = ['Telc(1)', 'Telc(2)', 'Tion(1)', 'Tion(2)']
    d = {}
    for tag in tags:
        d[tag] = gs.readDB(tag)
        for f in fields:
            plt.plot(d[tag]['Time']*scale, d[tag][f], '-', label="{}: {}".format(tag[5:],f))

    plt.xlabel('Time')
    plt.ylabel('Temp')
    if tmax > tmin: plt.xlim([tmin*scale, tmax*scale])
    if Tmax > Tmin: plt.ylim([Tmin, Tmax])
    plt.legend(loc=3, prop={'size':10}, bbox_to_anchor=(1,0)).draggable()
    plt.tight_layout(pad=7)
    plt.show()

def main(argv):
    """ Parse the command line and call the parsing and plotting routine """
    tmin  = 0.0
    tmax  = 0.0
    Tmin  = 0.0
    Tmax  = 0.0
    scale = 1.0
    
    helpString = '[--tmax=<tmax> | --Tmax <Tmax> | --tmin=<tmin> | --Tmin=<Tmin> | --scale=<scale> | -h | --help] <runtags>';

    try:
        opts, args = getopt.getopt(argv,"h", ["tmin=", "tmax=", "Tmin=", "Tmax=", "scale=", "help"])
    except getopt.GetoptError:
        print sys.argv[0], helpString
        sys.exit(2)
    for opt, arg in opts:
        if opt in ('-h', '--help'):
            print sys.argv[0], helpString
            sys.exit()
        elif opt in ("--tmin"):
            tmin = float(arg)
        elif opt in ("--tmax"):
            tmax = float(arg)
        elif opt in ("--Tmin"):
            Tmin = float(arg)
        elif opt in ("--Tmax"):
            Tmax = float(arg)
        elif opt in ("--scale"):
            scale = float(arg)

    if len(args)<=0:
        print sys.argv[0], helpString
        sys.exit(2)

    plot(args, tmin, tmax, Tmin, Tmax, scale)

if __name__ == "__main__":
    main(sys.argv[1:])
