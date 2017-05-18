
#!/usr/bin/python2

import sys, os, re
sys.path.append('')
os.environ['srcdir'] = '/home/herve/Documents/Script_GNURadio/IEEE_802.11/gr-ieee802-11/lib'
os.chdir('/home/herve/Documents/Script_GNURadio/IEEE_802.11/gr-ieee802-11/build/lib')

if __name__ == '__main__':
    sys.path.append('/home/herve/Documents/Script_GNURadio/IEEE_802.11/gr-ieee802-11/lib/../python')
    import build_utils
    root, inp = sys.argv[1:3]
    for sig in sys.argv[3:]:
        name = re.sub ('X+', sig, root)
        d = build_utils.standard_dict(name, sig, 'ieee802_11')
        build_utils.expand_template(d, inp, '_impl')
