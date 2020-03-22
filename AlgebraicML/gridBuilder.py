for i in range(16):
    bstr = "{0:b}".format(i).zfill(4)
    out = str(i) + ","
    for char in bstr:
        out += char + ","
    if (bstr[0] == "1" and bstr[2] == "1") or (bstr[1]=="1" and bstr[3]=="1"):
        out += "1"
    else:
        out += "0"
    print(out)