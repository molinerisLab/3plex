---
title: "R Notebook"
output: html_notebook
---

```{r}
z <- read.delim("pluto", header=T)
```


```{r}
peaks <- c("idr.conservative", "idr.optimal", "overlap.conservative.qTh005", "overlap.optimal.qTh005")

for (peak in peaks) {
  p <- ggboxplot(z, x = "technique", y = peak, color = "technique", palette = c("#00AFBB", "#E7B800", "#FC4E07")) + ggtitle(peak)
  print(p)
}
```

```{r}
for (peak in peaks) {
  p <- ggline(z, x = "technique", y = peak, add = c("mean_se", "jitter"))
  print(p)
}

```


```{r}
for (peak in peaks) {
  kt <- kruskal.test(z[,peak] ~ z[,"technique"])
  print(peak)
  print(kt)
}

```


```{r}
for (peak in peaks) {
  wt <- pairwise.wilcox.test(z[,peak], z[,"technique"], p.adjust.method = "BH")
  print(paste0("\n",peak,"\n"))
  print(wt)
}

```









