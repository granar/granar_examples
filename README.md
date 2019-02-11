# How to use GRANAR

GRANAR was built to reconstruct root cross-section. It can generate root from a simple set of parameters.

To run this example, the granar package is required. This package can be found on [GitHub](https://github.com/granar/granar).
Alternatively, by simply executing the following line in the R environment:

```{r}
install.packages("devtools")
library(devtools)
install_github("granar/grarar")
library(granar)
```

```{r preambule, echo=T, warning=F, message=F}
# Load the libraries
library(xml2)
library(tidyverse)
library(plyr)
library(deldir)
library(alphahull)
library(viridis)
library(cowplot)
library(retistruct)
library(Hmisc)
```
```{r, eval=T, echo= F, warning=FALSE, message=F}
library(granar)
```

```{r}
# Read the parameters
params <- read_param_xml("Zea_mays_Heymans_2019.xml")
```

```{r GRANAR, message = F, warning= F}
# Run the model
while (TRUE) {
  sim <- try(create_anatomy(parameters = params, verbatim=F), silent = TRUE)
  if(!is(sim,'try-error')) break
}
```

```{r plot}
# plot the output
plot_anatomy(sim, col = "type")
```
