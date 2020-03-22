for i in range(16):
    bstr = "{0:b}".format(i).zfill(4)
    out = str(i) + ","
    for char in bstr:
        out += char + ","
    if int(bstr[0])*int(bstr[3]) - int(bstr[1])*int(bstr[2]):
        out += "1"
    else:
        out += "0"
    print(out)