import csv
import numpy as np
import sys,os,re
import urllib2
import json






######################################################
def process_model_file(model_file):
    '''
    This function parses a csv model file into a dict to provide model simulation and return the dict and kmer lenght
    '''
    model_kmers = dict()
    with open(model_file, 'rb') as csv_file:
        reader = csv.reader(csv_file, delimiter="\t")
        d = list(reader)
        header = d[0]
        #print header
        meancolumnid = header.index("level_mean")
        kmercolumnid = header.index("kmer")

        d2 = np.array(d[1:])
        #for x in d2:
        #    print x
        model_kmers=dict(d2[:,[kmercolumnid,meancolumnid]])
        kmer_length=len(d2[0][kmercolumnid])
    return model_kmers,kmer_length



######################################################
def _urlopen(url, *args):

    """Open a URL, without using a proxy for localhost.

    While the no_proxy environment variable or the Windows "Bypass proxy

    server for local addresses" option should be set in a normal proxy

    configuration, the latter does not affect requests by IP address. This

    is apparently "by design" (http://support.microsoft.com/kb/262981).

    This method wraps urllib2.urlopen and disables any set proxy for

    localhost addresses.

    """
    #print "_urlopen called"
    try:

        host = url.get_host().split(':')[0]
    #    print "Yo mamma",host

    except AttributeError:

        host = urlparse.urlparse(url).netloc.split(':')[0]

    import socket

    #print "socket in"

    # NB: gethostbyname only supports IPv4

    # this works even if host is already an IP address

    addr = socket.gethostbyname(host)

    #if addr.startswith('127.'):
    #    print "returning no proxt"
    #    return _no_proxy_opener.open(url, *args)

    #else:
    #    print "returning urlib"
    return urllib2.urlopen(url, *args)


def execute_command_as_string( data, host = None, port = None):
        '''
        This function is used for sending messages to and from minKNOW.
        '''
        host_name = host

        port_number = port
        #print data

        url = 'http://%s:%s%s' % (host_name, port_number, "/jsonrpc")
        #print url

        req = urllib2.Request(url, data=data, headers={'Content-Length': str(len(data)), 'Content-Type': 'application/json'})
        f = None
        try:
            f = _urlopen(req)
        except Exception, err:
            err_string = "Failed to connect to minKNOW. Likely reasons include minKNOW not running, the wrong IP address for the minKNOW server or firewall issues. Two way control will not be possible with minKNOW."
            print err_string, err
        json_respond = json.loads(f.read())

        f.close()

        return json_respond

def send_message(message,ipadd):
    '''
    This function is used to send messages to and from minKNOW
    '''
    message_to_send = '{"id":"1", "method":"user_message","params":{"content":"%s"}}' %(message)
    results = execute_command_as_string(message_to_send,ipadd,8000)
    return results


def query_yes_no(question, default="yes"):
    """Ask a yes/no question via raw_input() and return their answer.

    "question" is a string that is presented to the user.
    "default" is the presumed answer if the user just hits <Enter>.
        It must be "yes" (the default), "no" or None (meaning
        an answer is required of the user).

    The "answer" return value is True for "yes" or False for "no".
    """
    valid = {"yes": True, "y": True, "ye": True,
             "no": False, "n": False}
    if default is None:
        prompt = " [y/n] "
    elif default == "yes":
        prompt = " [Y/n] "
    elif default == "no":
        prompt = " [y/N] "
    else:
        raise ValueError("invalid default answer: '%s'" % default)

    while True:
        sys.stdout.write(question + prompt)
        choice = raw_input().lower()
        if default is not None and choice == '':
            return valid[default]
        elif choice in valid:
            return valid[choice]
        else:
            sys.stdout.write("Please respond with 'yes' or 'no' "
                             "(or 'y' or 'n').\n")
