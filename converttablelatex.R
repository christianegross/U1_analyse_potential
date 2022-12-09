args <- commandArgs(trailingOnly = TRUE)
name <- as.character(args[1])
namenew <- as.character(args[2])

data <- read.table(name, sep = ",", header = TRUE)
print(data)
write.table(data, namenew, sep = "$ & $",
        eol = "$\\\\\n$", col.names = TRUE, row.names = FALSE)
string <- sprintf("the table %s has been copied to %s with latex style",
                    name, namenew)
print(string)