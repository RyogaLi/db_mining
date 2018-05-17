import requests
import socket


if __name__ == '__main__':
	session = requests.Session()
	session.proxies = {}
	session.proxies['http'] = 'socks5h://localhost:9050'
	session.proxies['https'] = 'socks5h://localhost:9050'
	r = session.get("http://www.google.com")
	print(r.text)