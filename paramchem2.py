from bs4 import BeautifulSoup
from xml.dom import minidom
import mechanize
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-u", "--username")
parser.add_argument("-p", "--password")
parser.add_argument("-c", "--conf")
args = parser.parse_args()

br = mechanize.Browser()
br.set_handle_robots(False)
br.set_handle_refresh(False)
br.addheaders = [('User-agent', 'Firefox')]
br.set_handle_redirect(True)

url = "https://cgenff.paramchem.org/userAccount/userLogin.php"
response = br.open(url)

br.form = list(br.forms())[0]
br["usrName"] = args.username
br["curPwd"] = args.password
response = br.submit()

filename = args.conf
br.form = list(br.forms())[0]
br.form.add_file(open(filename, "rb"), 'text/plain', filename)

response = br.submit()
xml = response.read().decode("utf-8").strip()

print(xml)

dom = minidom.parseString(xml)