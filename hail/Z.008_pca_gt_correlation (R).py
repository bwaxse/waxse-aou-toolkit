#!/usr/bin/env python
# coding: utf-8

# In[ ]:


# Install if needed
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("SNPRelate", "gdsfmt"))
install.packages(c("data.table", "ggplot2", "dplyr"))


# In[ ]:


library(SNPRelate)
library(data.table)
library(ggplot2)
library(dplyr)


# In[ ]:




