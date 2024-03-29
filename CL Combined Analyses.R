#Install packages
library(tidyverse)
library(agricolae)
library(car)
library(survival)
library(RColorBrewer)
library(scales)
library(viridis)

#Bioassay analysis
CL.dat <- read.csv("CL EHS bioassay data.csv", fileEncoding = 'UTF-8-BOM')

CL.dat2 <-  CL.dat %>%
  select(Day, Treatment, Replicate, Alive, Dead) %>%
  pivot_longer(cols=c(-Day, -Treatment, -Replicate), names_to="Status", values_to = "Individuals") %>%
  mutate(Status=as.factor(Status)) %>%
  group_by(Treatment, Replicate, Day) %>%
  mutate(Total=sum(Individuals)) %>%
  filter(Status=="Dead" & !is.na(Individuals)) %>%
  group_by(Replicate, Treatment) %>%
  mutate(Dead = Individuals - lag(Individuals, default = first(Individuals), order_by = Day)) %>%
  mutate(Dead=ifelse(Dead>0, Dead, 0)) %>%
  mutate(Percent=1-Individuals/Total)

CL.dat3 = CL.dat2 %>%
  group_by(Treatment, Day) %>%
  summarize(mn.ind=mean(Individuals),
            se.ind=sd(Individuals)/sqrt(length(Individuals)),
            mn.total=mean(Total),
            se.total=sd(Total)/sqrt(length(Total)),
            mn.per=mean(Percent),
            se.per=sd(Percent)/sqrt(length(Percent)),
            n=length(Percent))

palette1 <- c("#d67e1a", "#1a8ed6")

theme = theme_bw()+theme(text = element_text(size=25), axis.title.x = element_text(size=30), axis.text.x = element_text(size=20), axis.text.y = element_text(size=25), title = element_text(size=35), legend.title = element_text(size=25), legend.text = element_text(size=20), strip.background=element_rect(fill="white"))

limits.per=aes(ymin=mn.per-se.per, ymax=mn.per+se.per)

plt1=ggplot(CL.dat3, aes(x=Day, y=mn.per, color=Treatment))+
  scale_color_manual(values=palette1)+theme+geom_errorbar(limits.per)+
  geom_line(size=3)+ylab("Mean survival")+scale_y_continuous(label=percent, limits=c(0,1))+scale_x_continuous(breaks=0:13, minor_breaks=F)+
  theme(legend.position=c(0.7,0.3))+theme(legend.text = element_text(size=30), axis.text.x = element_text(size=35), axis.text.y=element_text(size=35), axis.title.x=element_text(size=40), axis.title.y=element_text(size=40),
                                          panel.grid.major = element_blank(), panel.grid.minor = element_blank())

plt1

CL.dat4 = CL.dat3 %>%
  select(Treatment, Day, mn.ind, mn.total) %>%
  group_by(Treatment, Day) %>%
  pivot_longer(cols=c(mn.ind, mn.total), names_to="Status", values_to="Mean") %>%
  mutate(Status=replace(Status, Status=="mn.ind", "Dead"),
         Status=replace(Status, Status=="mn.total", "Alive"))

palette2 <- c("#d6c91a", "#3d3a08")

plt2=ggplot(CL.dat4, aes(x=Day, y=Mean, fill=Status))+
  scale_fill_manual(values=palette2)+ theme+
  geom_area(position="dodge")+facet_wrap(~Treatment)+ylab("Mean number of crawlers")+
  scale_x_continuous(breaks=0:13, minor_breaks=F)+
  theme(legend.position=c(0.1, 0.9), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.ticks=element_line(size=2))

plt2

CL.dat.mort= CL.dat2 %>%
  select(Day, Treatment, Dead, Replicate, Total) %>%
  rename(Mortality = Dead)

CL.dat.mort2 = CL.dat.mort %>%
  filter(Day==13) %>%
  mutate(Mortality=Total) %>%
  mutate(Day=14)

CL.dat.mort3=bind_rows(CL.dat.mort, CL.dat.mort2)

convert_mortality=function(mort.dat){
  dat=data.frame()
  for(i in 1:length(mort.dat$Mortality)){
    Day=mort.dat[i,]$Day
    Treatment=mort.dat[i,]$Treatment
    Replicate=mort.dat[i,]$Replicate
    Mortality=mort.dat[i,]$Mortality
    if(Mortality>0){
      temp=data.frame(Day=rep(Day, times=Mortality),
                      Treatment=rep(Treatment, times=Mortality),
                      Replicate=rep(Replicate, times=Mortality),
                      Status=rep(ifelse(Day!=14, 1, 0), times=Mortality))
      dat=bind_rows(dat, temp)
    }
  }
  return(dat)
}

KM.dat=convert_mortality(CL.dat.mort3)

survdiff(Surv(KM.dat$Day, KM.dat$Status)~KM.dat$Treatment)

#Morphogroup analysis
n.dat <- read.csv("Needle sterilization groups.csv")

n.dat2 <- n.dat %>%
  mutate(Needle.status=as.factor(Needle.status), EHS.Status=as.factor(EHS.Status), Fungus.status=as.factor(Fungus.status), Bleach=as.factor(Bleach), Site=as.factor(Site)) %>%
  mutate(EHS.Status=recode_factor(EHS.Status, Negative="EHS-", Positive="EHS+")) %>%
  mutate(Fungus.status=recode_factor(Fungus.status, Negative="CL-", Positive="CL+")) %>%
  mutate(Bleach=recode_factor(Bleach, No="Unsterilized", Yes="")) %>%
  mutate(Needle.status=recode_factor(Needle.status, Dead=0.6, Alive=1)) %>%
  mutate(Site=recode_factor(Site, BAF="Deal Family Farm", MR="Mount Rogers", MRO="Mount Rogers Orchard",
                            UM="Upper Mountain", VF="Vannoy Farm")) %>%
  group_by(Needle.status, EHS.Status, Fungus.status, Bleach, Site, Group.Number) %>%
  tally()

pal=colorRampPalette(brewer.pal(8, "Dark2"))(31)

plt3=ggplot(n.dat2, aes(x=paste(EHS.Status, Fungus.status, Bleach, sep="\n"), y=n, alpha=as.numeric(Needle.status), fill=as.factor(Group.Number)))+
  geom_col(position=position_dodge2(preserve = "single"))+facet_wrap(~Site)+scale_fill_manual(values=pal)+
  scale_y_continuous(breaks=seq(0, 5, 1))+theme+theme(legend.position=c(0.83,0.2))+
  theme(axis.title.x=element_blank())+ylab("Count")+labs(fill="Morphogroup")+guides(fill=guide_legend(ncol=4), alpha="none")

plt3

#Spore measurement analysis
spo <- read.csv("CL isolate measurements.csv", fileEncoding = 'UTF-8-BOM')

spo2 = spo %>%
  group_by(Strain, Dye, Measurement) %>%
  summarize(mn=mean(Value), sd=sd(Value), min=min(Value), max=max(Value), n=length(Value))

spo_comb = spo %>%
  group_by(Dye, Measurement) %>%
  summarize(mn=mean(Value), sd=sd(Value), min=min(Value), max=max(Value), n=length(Value)) %>%
  mutate(Strain="Combined")

sp.mod=aov(Value~Strain+Measurement, data=spo)
summary(sp.mod)

shapiro.test(rstandard(sp.mod))

spo.width=spo %>% filter(Measurement=="width")
with(spo.width, leveneTest(Value, factor(Strain)))

spo.length=spo %>% filter(Measurement=="length")
with(spo.length, leveneTest(Value, factor(Strain)))?

sp.tuk=HSD.test(sp.mod, trt=c("Strain", "Measurement"))

sp.tuk.r=sp.tuk$means %>%
  rownames_to_column() %>%
  separate(rowname, into=c("Strain", "Measurement"), sep=":") %>%
  select(Strain, Measurement, r)

sp.grouping=sp.tuk$groups %>%
  rownames_to_column() %>%
  separate(rowname, into=c("Strain", "Measurement"), sep=":") %>%
  left_join(sp.tuk.r) %>%
  rename("n"=r)

sp.grouping

#write.csv(sp.grouping, "spore groupings.csv", row.names=F)

spo2_comb = spo2 %>%
  bind_rows(spo_comb) %>%
  mutate(Strain=as.factor(Strain)) %>%
  mutate(Measurement=str_replace(Measurement, "length", "Length"),
         Measurement=str_replace(Measurement, "width", "Width"),)

spo2_dyes = spo2 %>%
  filter(Strain=="VF 11") %>%
  mutate(Strain=as.factor(Strain)) %>%
  mutate(Measurement=str_replace(Measurement, "length", "Length"),
         Measurement=str_replace(Measurement, "width", "Width"),)

spo3 = spo2 %>%
  pivot_longer(cols=c(mn, sd, max, min, n), values_to="Value", names_to="Stat") %>%
  pivot_wider(id_cols=c(Strain, Dye), names_from=c(Measurement, Stat), values_from=Value)


#Temperature growth analysis
growth_dat <- read.csv("CL Temperature Assay.csv")

growth2 <-  growth_dat %>%
  select(Site, Site.Isolate, Replicate, Timepoint, Top, Bottom, Left, Right) %>%
  pivot_longer(cols=c(Top, Bottom, Left, Right), names_to="Direction", values_to = "Radius") %>%
  group_by(Site, Site.Isolate, Replicate, Timepoint) %>%
  summarize(avg.Rad=mean(Radius))

growth_stat <- growth2 %>%
  group_by(Timepoint, Site, Site.Isolate) %>%
  summarise(Culture_radius=mean(avg.Rad),
            n=length(avg.Rad),
            SE=sd(avg.Rad)/sqrt(n)
  ) %>%
  mutate(Site.num=str_c(Site, Site.Isolate))

growth_rate= growth_stat %>%
  ungroup() %>%
  filter(Site!="UM") %>%
  group_by(Site, Site.Isolate) %>%
  mutate(Rate = Culture_radius - lag(Culture_radius, default = first(Culture_radius))) %>%
  filter(Timepoint>2) %>%
  ungroup() %>%
  summarize(avg_rate=mean(Rate)/2)

growth_rate

#Published and observed spore metrics

con_dat <- read.csv("Conoideocrella spore measurements.csv")

con_dat2 = con_dat %>%
  rename_with(str_to_title) %>%
  mutate(Metric=str_to_title(Metric)) %>%
  pivot_longer(cols=Mean:Max, names_to="Metric2", values_to="Value") %>%
  pivot_wider(names_from=c("Metric", "Metric2"), values_from="Value") %>%
  mutate(Species=as.factor(Species))

CL_limits=con_dat2 %>%
  mutate(CL=ifelse(grepl("is", Species), "No", "Yes")) %>%
  group_by(CL) %>%
  summarize(Length_Max=max(Length_Mean+Length_Difference, na.rm=T), Length_Min=min(Length_Mean-Length_Difference, na.rm=T), Width_Max=max(Width_Mean+Width_Difference, na.rm=T), Width_Min=min(Width_Mean-Width_Difference, na.rm=T)) %>%
  filter(CL=="Yes") %>%
  select(-CL)

theme = theme_bw()+theme(text = element_text(size=15), axis.title.x = element_text(size=20), axis.title.y = element_text(size=20), axis.text.x = element_text(size=20), axis.text.y = element_text(size=20), title = element_text(size=35), legend.title = element_text(size=15))

plt4=ggplot(con_dat2)+
  geom_rect(aes(ymin=Length_Mean-Length_Difference, ymax=Length_Mean+Length_Difference, xmin=Width_Mean-Width_Difference, xmax=Width_Mean+Width_Difference, fill=Species), alpha=0.5)+
  geom_point(aes(x=Width_Mean, y=Length_Mean, color=Species), size=4)+
  geom_point(aes(x=Width_Mean, y=Length_Mean), color="white", size=2)+
  theme_bw()+scale_color_viridis_d(option="C")+
  scale_fill_viridis_d(option="C")+theme+ylim(4,16)+xlim(1,5)+
  theme(text=element_text(size=15))+ylab("Length")+xlab("Width")+
  ggrepel::geom_label_repel(aes(x=Width_Mean, y=Length_Mean, label=Species), color="black", box.padding=1)+
  guides(color="none")+
  theme(axis.ticks.length=unit(-0.25, "cm"),
        legend.position="bottom", panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  geom_rect(dat=CL_limits, aes(ymin=Length_Min, ymax=Length_Max,
                               xmin=Width_Min, xmax=Width_Max),
            color="black",fill=NA, linewidth=1, linetype=2)+
  ggrepel::geom_label_repel(dat=CL_limits, aes(y=Length_Min-0.3, x=Width_Max-0.1),
            label="Conoideocrella luteorostrata")

plt4
