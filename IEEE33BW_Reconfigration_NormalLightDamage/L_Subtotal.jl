#Subtotal v1.0
using DataStructures 
for term in Tuple(readdir(".//report//"))
    listFile=Tuple(readdir(".//report//" * term))
    if (term*"_summary.txt" in listFile)==false
        numPoint=[]
        for fileName in listFile
            temp=split(fileName)
            push!(numPoint,parse(Int64,temp[3]))
        end
        fileDic=SortedDict(Tuple(numPoint).=>listFile)
        fileTerm=open("report//"*term*"//"*term*"_summary.txt","a+")
        for num in keys(fileDic)
            write(fileTerm,readlines(".//report//"*term*"//"*fileDic[num], keep=true)[1])
        end
        close(fileTerm)
    end
end