files = readdir()
deletefiles = filter(files) do f
endswith(f, ".png") || endswith(f, ".bin") || endswith(f, ".txt")
end
foreach(deletefiles) do f
    if isfile(f)
        rm(f)
    end
end
