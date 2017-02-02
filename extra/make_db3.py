import requests

def get_spectra_from_file():
    f = open('metadata')
    spectra_dict = {}
    for line in f.readlines():
        d = eval(line.strip())
        spectra_dict[d['SpectrumID']] = d

    f.close()
    urls = []
    for spectrum in spectra_dict:
        SERVER_URL = "http://gnps.ucsd.edu/ProteoSAFe/SpectrumCommentServlet?SpectrumID="
        url = SERVER_URL + spectrum
        urls.append(url)
    return urls

print "reading metadata"
urls = get_spectra_from_file()

f = open('all_data', 'a+')

i = 0
err = 0
for url in urls:
    try:
        i += 1
        print i
        data1 = requests.get(url).text
        f.write(data1 + "\n")
    except:
        err += 1
        print "ERRRR"

f.close()
print "Number of errors:"
print err