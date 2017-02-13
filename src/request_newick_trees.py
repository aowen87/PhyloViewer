#!/usr/bin/python
'''
Submit requests to phyloT to generate
a bunch of phylogenetic trees in the form
of newick strings. 

'''
import mechanize
import cookielib
import argparse
import glob


if __name__ == "__main__":

    br = mechanize.Browser()
    parser = argparse.ArgumentParser()
    parser.add_argument('input_dir')
    parser.add_argument('out_path')
    args     = parser.parse_args()
    out_path = args.out_path
    in_path  = args.input_dir
    
    # Cookie Jar
    cj = cookielib.LWPCookieJar()
    br.set_cookiejar(cj)

    # Browser options
    br.set_handle_equiv(True)
    br.set_handle_redirect(True)
    br.set_handle_referer(True)
    br.set_handle_robots(False)
    br.set_handle_refresh(mechanize._http.HTTPRefreshProcessor(), max_time=1)

    if out_path[-1] != '/':
        out_path = out_path + '/'

    if in_path[-1] != '/':
        in_path = in_path + '/'

    in_grab = in_path + '*'

    for pth in glob.glob(in_grab):
        #I might not have to open this page
        #for every run, but I'm not sure how 
        #to avoid doing that...  
        br.open("http://phylot.biobyte.de/")
        br.select_form(nr=0)
        sample_num = pth.split('/')[-1].split('_')[-1]
        in_f   = open(pth, 'r')
        sample = in_f.readline()
        br.form['binary'] = ['1']
        br.form['treeElements'] = sample
        response = br.submit()
        out_name = out_path + "tree_" + sample_num
        out_f = open(out_name, 'w+')
        out_f.write(response.read())
        out_f.close()



