library(tidyr)
library(dplyr)

# 1. Load your tidy/long data
df_long <- read.csv("TIDY_Transpose of ImmuneGenes.csv")

# 2. Pivot to wide format
# names_from: the column that contains gene names (to become headers)
# values_from: the column that contains the actual expression values
df_wide <- df_long %>%
  pivot_wider(names_from = Gene, values_from = fpkm)

# 3. Sort by Day (optional, but helpful for the NSV script)
df_wide <- df_wide %>% arrange(Day)

# 4. Save the wide file for use in the NSV script
write.csv(df_wide, "ImmuneGenes_WIDE.csv", row.names = FALSE)