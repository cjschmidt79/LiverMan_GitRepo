# 1. Load the libraries
if (!require("tidyr")) install.packages("tidyr")
if (!require("dplyr")) install.packages("dplyr")

library(tidyr)
library(dplyr)

# 2. Read the file you already have
# Ensure "Transcriptome_Tidy_Long4columnsFORNSV.csv" is in your working directory
df_long <- read.csv("Transcriptome_Tidy_Long4columnsFORNSV.csv")

# 3. Pivot from Long to Wide
# This moves the 'Gene' names into columns and 'fpkm' into the values
df_wide <- df_long %>%
  pivot_wider(names_from = Gene, values_from = fpkm)

# 4. Save the new file
write.csv(df_wide, "WIDE_LiverMetabolome.csv", row.names = FALSE)

cat("Success! 'WIDE_LiverMetabolome.csv' has been created.\n")