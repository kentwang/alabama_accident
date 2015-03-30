library(ISLR)                          # to access Hitters data
player = rownames(Hitters)             # player info stored in row names
player = gsub('-','',player)           # remove the leading '-' in player name
H = data.frame(player,Hitters)         # add player to the data
H = subset(H, HmRun>10)                # only use players with at least 10 Homeruns
H$avg = H$Hits/H$AtBat 

#- Try with barplot (order from most to least)
ggplot(H,aes(x=reorder(player,HmRun, FUN = function(V) v),y=HmRun,fill=avg>=.300)) + geom_bar(stat="identity") +  coord_flip() + theme(axis.text.y=element_text(size=rel(0.8)))
#> Warning: position_stack requires constant width: output may be incorrect

#- Use geom_point to make the dot plot
dot_theme = theme_bw() +
  theme(panel.grid.major.x=element_blank(),
        panel.grid.minor.x=element_blank(),
        panel.grid.major.y=element_line(color="grey60",
                                        linetype="dashed"))

ggplot(H,aes(y=reorder(player, HmRun, FUN = function(x) x),x=HmRun,color=avg>=.300)) + 
  geom_point(aes(size=avg)) + 
  dot_theme +
  geom_point(H,aes(y=reorder(player, HmRun, FUN = function(x) x),x=HmRun,color=avg>=.300))


H1 = H
H1$HmRun = H1$HmRun - 10

ggplot(rbind(data.frame(H, group = "a"), data.frame(H1, group = "b")),
       aes(y=reorder(player, HmRun, FUN = function(x) x),x=HmRun,color=avg>=.300)) + 
  geom_point(aes(size=avg, shape = group)) + 
  dot_theme