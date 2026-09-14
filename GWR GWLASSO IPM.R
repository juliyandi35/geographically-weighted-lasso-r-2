# Import peta
library(sf)
Ind_map = read_sf('Peta Kabupaten/BATAS KABUPATEN KOTA DESEMBER 2019 DUKCAPIL.shp')
Kab_Jatim <- c("PACITAN","PONOROGO","TRENGGALEK","TULUNGAGUNG",
               "BLITAR","KEDIRI","MALANG","LUMAJANG","JEMBER",
               "BANYUWANGI","BONDOWOSO","SITUBONDO","PROBOLINGGO",
               "PASURUAN","SIDOARJO","MOJOKERTO","JOMBANG","NGANJUK",
               "MADIUN","MAGETAN","NGAWI","BOJONEGORO","TUBAN",
               "LAMONGAN","GRESIK","BANGKALAN","SAMPANG","PAMEKASAN",
               "SUMENEP","KOTA KEDIRI","KOTA BLITAR","KOTA MALANG","KOTA PROBOLINGGO",
               "KOTA PASURUAN","KOTA MOJOKERTO","KOTA MADIUN","KOTA SURABAYA","KOTA BATU")
Jatim_map <- subset(Ind_map, KAB_KOTA %in% Kab_Jatim)
Jatim_map

# Membuat plot peta
library(ggplot2)
ggplot() +
  geom_sf(data = Jatim_map) +
  labs(title = "Peta Kabupaten/Kota di Jawa Timur")

#Statistika Deskriptif
library(readxl)
IPM1<-read_excel("IPM GWR.xlsx")
colnames(IPM1) <- c("KAB_KOTA","LONG","LAT","Y","X1","X2","X3","X4")
summary(IPM1)

Data <- merge(Jatim_map,IPM1,by="KAB_KOTA")

library(sp)
Data <- st_as_sf(Data)
ggplot(Data) + geom_sf(aes(fill = Y)) +
  scale_fill_gradient2(
    midpoint = 73, low = "#28E2E5", high = "#DF536B") +
  theme_bw()+
  geom_text(
    aes(label = KAB_KOTA, x = coordinates(as(Data,"Spatial"))[,1], y = coordinates(as(Data,"Spatial"))[,2]),
    vjust = -0.5,
    color = "black",
    size = 1.5,
    check_overlap = TRUE
  )+ggtitle("IPM di Jawa Timur")+xlab("Longitude")+ylab("Latitude")

# Model OLS
regresi<-lm(formula=Y~X1+X2+X3+X4,data=IPM1)
summary(regresi)

#Uji Heterogenitas (Breusch-Pagan)
library(lmtest)
bptest(regresi)

# Uji Multikolinearitas
library(car)
vif(regresi)

library(GWmodel)
library(spdep)
#Penskalaan Data ke Dalam Data Frame Spasial# 
library(dplyr)

# moran test
x = coordinates(as(Data,"Spatial"))[,1]
y = coordinates(as(Data,"Spatial"))[,2]
coords <-cbind(x,y)
jarak <-as.matrix(1/dist(coords))
lm.morantest(regresi,listw=mat2listw(jarak), alternative="two.sided") # Moran Index

# GWR
# convert to sp
Data.sp = as(Data, "Spatial")
bw<-bw.gwr(Y~X1+X2+X3+X4, data=Data.sp,
                 approach="CV",kernel="exponential",adaptive=FALSE,
                 p=2,longlat=FALSE)

#Matriks Pembobot model GWR#
library(fields)
W<-exp(-jarak/bw)
W
#Pemodelan GWR 
GWRModel<-gwr.basic(Y~ X1+X2+X3+X4,data=Data.sp,
                          bw=bw,kernel="exponential",adaptive=FALSE,p=2,
                          F123.test=T)
# Uji F
GWRModel$Ftests$F1.test
GWRModel$Ftests$F2.test
GWRModel$Ftests$F3.test
GWRModel$Ftests$F4.test

#Estimasi parameter model GWR
GWR_sf = st_as_sf(GWRModel$SDF)

# T hitung dengan alfa = 0.05 dan DF = 38-5 adalah 2.034515 
# Tentukan daerah signifikan
Data$tval.X1 = GWR_sf$X1_TV
Data$pvalue.X1 = 2*pt(abs(Data$tval.X1), df = 2.034515, lower.tail = FALSE)

Data$tval.X2 = GWR_sf$X2_TV
Data$pvalue.X2 = 2*pt(abs(Data$tval.X2), df = 2.034515, lower.tail = FALSE)

Data$tval.X3 = GWR_sf$X3_TV
Data$pvalue.X3 = 2*pt(abs(Data$tval.X3), df = 2.034515, lower.tail = FALSE)

Data$tval.X4 = GWR_sf$X4_TV
Data$pvalue.X4 = 2*pt(abs(Data$tval.X4), df = 2.034515, lower.tail = FALSE)

# Pakai file excel yang beda
df.lasso = read_excel("IPM Skripsi.xlsx")
colnames(df.lasso) <- c("KAB_KOTA","LONG","LAT","Y","X1","X2","X3","X4","X5","X6","X7")
Data.lasso <- merge(Jatim_map,df.lasso,by="KAB_KOTA")
df <- st_drop_geometry(Data.lasso)
names(df)
df <- df[,-c(1:3)]
Check.NA <- colSums(is.na(df))
Check.NA

# Model OLS
OLS.LASSO<-lm(formula=Y~X1+X2+X3+X4+X5+X6+X7,data=df.lasso)
summary(OLS.LASSO)

#Uji Heterogenitas (Breusch-Pagan)
library(lmtest)
bptest(OLS.LASSO)

# Uji Multikolinearitas
library(car)
vif(OLS.LASSO)

## GWLASSO
split_value = 0.8
# Step2: Split data into training and testing sets
set.seed(22)
train_idx <- sample(nrow(df), split_value * nrow(df))
train_data <- df[train_idx,]
nrow(train_data)
test_data <- df[-train_idx, ]
nrow(test_data)

names(train_data)

x_train <-train_data[,-1]
y_train <-train_data[,1]
x_test <- test_data[,-1]
y_test <- test_data[,1]

# Step3: Important variable selection using LASSO
# fit LASSO model for variable selection

train_data_LASSO<-cbind(y_train,x_train)
x_lasso <- model.matrix(y_train ~ ., train_data_LASSO)
y_lasso <- train_data_LASSO[, "y_train"]
y_lasso <- y_lasso[1:nrow(x_lasso)]

### select important variables based on LASSO model
library(glmnet)
set.seed(123)
cvfit <- cv.glmnet(x_lasso, y_lasso, alpha = 1, standardize = TRUE);plot(cvfit)
coef_path <- predict(cvfit, type="coefficients", s=cvfit$lambda.min)
selected_vars <- row.names(coef_path)[which(coef_path != 0)[-1]] # exclude intercept
coef_path.df <- data.frame(Variables = rownames(coef_path),Values = coef_path[,1])
coef_path.df <- coef_path.df[which(coef_path.df$Values!=0),]
coef_path.df <- coef_path.df[-1,] # exclude intercept
coef_path.df

OLS_Post.LASSO <-lm(formula=Y~X1+X3+X6,data=df.lasso)
summary(OLS_Post.LASSO)

#Uji Heterogenitas (Breusch-Pagan)
library(lmtest)
bptest(OLS_Post.LASSO)

# Uji Multikolinearitas
library(car)
vif(OLS_Post.LASSO)

# Extracting important variables
imp_var_lasso <- x_train[, selected_vars]
names(imp_var_lasso)

# extract optimal lambda value
opt_lambda <- cvfit$lambda.min
opt_lambda

# GWLASSO
# convert to sp
Data.lasso.sp = as(Data.lasso, "Spatial")
bw.lasso<-bw.gwr(Y~X1+X3+X6, data=Data.lasso.sp,
           approach="CV",kernel="exponential",adaptive=FALSE,
           p=2,longlat=FALSE)

#Pemodelan GWR 
GWLASSOModel<-gwr.basic(Y~ X1+X3+X6,data=Data.lasso.sp,
                    bw=bw,kernel="exponential",adaptive=FALSE,p=2,
                    F123.test=T)
# Uji F
GWLASSOModel$Ftests$F1.test
GWLASSOModel$Ftests$F2.test
GWLASSOModel$Ftests$F3.test
GWLASSOModel$Ftests$F4.test

#Estimasi parameter model GWR
GWLASSO_sf = st_as_sf(GWLASSOModel$SDF)

# T hitung dengan alfa = 0.05 dan DF = 38-3 adalah 2.032244509

# Tentukan daerah signifikan
Data.lasso$tval.X1 = GWLASSO_sf$X1_TV
Data.lasso$pvalue.X1 = 2*pt(abs(Data.lasso$tval.X1), df = 2.032244509, lower.tail = FALSE)

Data.lasso$tval.X3 = GWLASSO_sf$X3_TV
Data.lasso$pvalue.X3 = 2*pt(abs(Data.lasso$tval.X3), df = 2.032244509, lower.tail = FALSE)

Data.lasso$tval.X6 = GWLASSO_sf$X6_TV
Data.lasso$pvalue.X6 = 2*pt(abs(Data.lasso$tval.X6), df = 2.032244509, lower.tail = FALSE)

#Pemilihan Model 
# RMSE
RMSE <- function(actual,predicted){
  sqrt(mean((actual-predicted)^2))
}
OLS_RMSE <- RMSE(Data$Y,predict(regresi))
GWR_RMSE <- RMSE(Data$Y,GWR_sf$yhat)
GWLASSO_RMSE <- RMSE(Data$Y,GWLASSO_sf$yhat)

Diagnostic_comparison_tabel <-
  cbind(c("OLS","GWR","GWLASSO"), 
        rbind(summary(regresi)$r.squared,GWRModel$GW.diagnostic$gw.R2,GWLASSOModel$GW.diagnostic$gw.R2),
        rbind(OLS_RMSE,GWR_RMSE,GWLASSO_RMSE))
colnames(Diagnostic_comparison_tabel) <- c("Metode","R Squared","RMSE")
rownames(Diagnostic_comparison_tabel) <- NULL
Diagnostic_comparison_tabel <- data.frame(Diagnostic_comparison_tabel)
Diagnostic_comparison_tabel

#---------------------------------------------------------------#
#    signfikansi variabel (variabel X1)
#---------------------------------------------------------------#
Data$signfikansi_X1 <- NA
# Signifikan
Data[Data$pvalue.X1 <= 0.05, "signfikansi_X1"] <- "Signifikan"

# Tidak Signifikan
Data[Data$pvalue.X1 > 0.05, "signfikansi_X1"] <- "Tidak Signifikan"

#------------------------------------------------
#Gabung data GWR dengan SHP
ggplot(data=Data) +
  geom_sf(mapping=aes(fill =signfikansi_X1)) +
  scale_fill_manual(values = c("#28E2E5", "#DF536B"))+
  labs(fill="signfikansi")+
  geom_text(
    aes(label = KAB_KOTA, x = coordinates(as(Data.lasso,"Spatial"))[,1], y = coordinates(as(Data.lasso,"Spatial"))[,2]),
    vjust = -0.5,
    color = "black",
    size = 1.5,
    check_overlap = TRUE
  )+ggtitle("Signifikansi X1")

#---------------------------------------------------------------#
#    signfikansi variabel (variabel X2)
#---------------------------------------------------------------#
Data$signfikansi_X2 <- NA
# Signifikan
Data[Data$pvalue.X2 <= 0.05, "signfikansi_X2"] <- "Signifikan"

# Tidak Signifikan
Data[Data$pvalue.X2 > 0.05, "signfikansi_X2"] <- "Tidak Signifikan"

#------------------------------------------------
#Gabung data GWR dengan SHP
ggplot(data=Data) +
  geom_sf(mapping=aes(fill =signfikansi_X2)) +
  scale_fill_manual(values = c("#DF536B"))+
  labs(fill="signfikansi")+
  geom_text(
    aes(label = KAB_KOTA, x = coordinates(as(Data.lasso,"Spatial"))[,1], y = coordinates(as(Data.lasso,"Spatial"))[,2]),
    vjust = -0.5,
    color = "black",
    size = 1.5,
    check_overlap = TRUE
  )+ggtitle("Signifikansi X2")

#---------------------------------------------------------------#
#    signfikansi variabel (variabel X3)
#---------------------------------------------------------------#
Data$signfikansi_X3 <- NA
# Signifikan
Data[Data$pvalue.X3 <= 0.05, "signfikansi_X3"] <- "Signifikan"

# Tidak Signifikan
Data[Data$pvalue.X3 > 0.05, "signfikansi_X3"] <- "Tidak Signifikan"

#------------------------------------------------
#Gabung data GWR dengan SHP
ggplot(data=Data) +
  geom_sf(mapping=aes(fill =signfikansi_X3)) +
  scale_fill_manual(values = c("#DF536B"))+
  labs(fill="signfikansi")+
  geom_text(
    aes(label = KAB_KOTA, x = coordinates(as(Data.lasso,"Spatial"))[,1], y = coordinates(as(Data.lasso,"Spatial"))[,2]),
    vjust = -0.5,
    color = "black",
    size = 1.5,
    check_overlap = TRUE
  )+ggtitle("Signifikansi X3")

#---------------------------------------------------------------#
#    signfikansi variabel (variabel X4)
#---------------------------------------------------------------#
Data$signfikansi_X4 <- NA
# Signifikan
Data[Data$pvalue.X4 <= 0.05, "signfikansi_X4"] <- "Signifikan"

# Tidak Signifikan
Data[Data$pvalue.X4 > 0.05, "signfikansi_X4"] <- "Tidak Signifikan"

#------------------------------------------------
#Gabung data GWR dengan SHP
ggplot(data=Data) +
  geom_sf(mapping=aes(fill =signfikansi_X4)) +
  scale_fill_manual(values = c("#28E2E5", "#DF536B"))+
  labs(fill="signfikansi")+
  geom_text(
    aes(label = KAB_KOTA, x = coordinates(as(Data.lasso,"Spatial"))[,1], y = coordinates(as(Data.lasso,"Spatial"))[,2]),
    vjust = -0.5,
    color = "black",
    size = 1.5,
    check_overlap = TRUE
  )+ggtitle("Signifikansi X4")

#---------------------------------------------------------------#
#    Plotting variabel signifikan
#---------------------------------------------------------------#
Significant_Map <- Data

# Buat kolom baru untuk kombinasi signfikansi
Significant_Map <- Significant_Map %>%
  mutate(Variabel_Signifikan = case_when(
    # Utuh
    signfikansi_X1 == "Signifikan" & signfikansi_X2 == "Signifikan" &
      signfikansi_X3 == "Signifikan" & signfikansi_X4 == "Signifikan"  ~ "X1,X2,X3,X4",
    
    # Eliminasi 1
    signfikansi_X1 == "Signifikan" & signfikansi_X2 == "Signifikan" &
      signfikansi_X3 == "Signifikan" ~ "X1,X2,X3",
    signfikansi_X1 == "Signifikan" & signfikansi_X2 == "Signifikan" &
      signfikansi_X4 == "Signifikan"  ~ "X1,X2,X4",
    signfikansi_X1 == "Signifikan" &
      signfikansi_X3 == "Signifikan" & signfikansi_X4 == "Signifikan"  ~ "X1,X3,X4",
    signfikansi_X2 == "Signifikan" &
      signfikansi_X3 == "Signifikan" & signfikansi_X4 == "Signifikan"  ~ "X2,X3,X4",
    # Eliminasi 2
    signfikansi_X1 == "Signifikan" & signfikansi_X2 == "Signifikan"  ~ "X1,X2",
    signfikansi_X1 == "Signifikan" & signfikansi_X4 == "Signifikan"  ~ "X1,X4",
    signfikansi_X3 == "Signifikan" & signfikansi_X4 == "Signifikan"  ~ "X3,X4",
    signfikansi_X2 == "Signifikan" & signfikansi_X3 == "Signifikan"  ~ "X2,X3",
    
    # Eliminasi 3
    signfikansi_X1 == "Signifikan"  ~ "X1",
    signfikansi_X2 == "Signifikan"  ~ "X2",
    signfikansi_X3 == "Signifikan"  ~ "X3",
    signfikansi_X4 == "Signifikan"  ~ "X4",
    
    TRUE ~ "Tidak Signifikan"
  ))

# Buat skema warna kustom
warna_custom <- c(
  "X1,X2,X3,X4" = "black",
  "X1,X2,X3" = "red",
  "X1,X2,X4" = "green",
  "X1,X3,X4" = "blue",
  "X2,X3,X4" = "cyan",
  "X1,X2" = "magenta",
  "X1,X4" = "yellow",
  "X3,X4" = "gray",
  "X2,X3" = "darkgray",
  "X1" = "lightgray",
  "X2" = "orange",
  "X3" = "brown",
  "X4" = "pink",
  "Tidak Signifikan" = "white"
)

ggplot(data = Significant_Map) +
  geom_sf(mapping=aes(geometry = geometry,fill = Variabel_Signifikan)) +
  scale_fill_manual(values = warna_custom)+
  labs(fill="Variabel Signifikan")+
  geom_text(
    aes(label = KAB_KOTA, x = coordinates(as(Data.lasso,"Spatial"))[,1], y = coordinates(as(Data.lasso,"Spatial"))[,2]),
    vjust = -0.5,
    color = "black",
    size = 1.5,
    check_overlap = TRUE
  )+ggtitle(("Variabel Signifikan tiap Lokasi"))+xlab("Longitude")+ylab("Latitude")
