#Buatlah jarngan syarat untuk mengenali 7 huruf tersebut masing-masing dengan menggunakan : 
#a.	Jaringan HEBB
#b.	Jaringan Perceptron
#c.	Jaringan Propagasi balik error
#d.	SOM 
#e.	LVQ
#menentukan variabel precditor
a1<-c(0,0,1,1,0,0,0,0,
      0,0,0,1,0,0,0,0,
      0,0,0,1,0,0,0,0,
      0,0,1,0,1,0,0,0,
      0,0,1,0,1,0,0,0,
      0,0,1,0,1,0,0,0,
      0,1,1,1,1,1,0,0,
      0,1,0,0,0,1,0,0,
      1,1,1,0,1,1,1,0)
b1<-c(1,1,1,1,1,1,1,0,
      0,1,0,0,0,0,0,1,
      0,1,0,0,0,0,0,1,
      0,1,0,0,0,0,0,1,
      0,1,1,1,1,1,1,0,
      0,1,0,0,0,0,0,1,
      0,1,0,0,0,0,0,1,
      0,1,0,0,0,0,0,1,
      1,1,1,1,1,1,1,0)

c1<-c(0,0,1,1,1,1,1,0,
      0,1,0,0,0,0,0,1,
      1,0,0,0,0,0,0,0,
      1,0,0,0,0,0,0,0,
      1,0,0,0,0,0,0,0,
      1,0,0,0,0,0,0,0,
      1,0,0,0,0,0,0,0,
      0,1,0,0,0,0,0,1,
      0,0,1,1,1,1,1,0)

d1<-c(1,1,1,1,1,1,1,0,
      0,1,0,0,0,0,0,1,
      0,1,0,0,0,0,0,1,
      0,1,0,0,0,0,0,1,
      0,1,0,0,0,0,0,1,
      0,1,0,0,0,0,0,1,
      0,1,0,0,0,0,0,1,
      0,1,0,0,0,0,0,1,
      1,1,1,1,1,1,1,0)

e1<-c(1,1,1,1,1,1,1,1,
     0,1,0,0,0,0,0,1,
     0,1,0,0,0,0,0,0,
     0,1,0,1,0,0,0,0,
     0,1,1,1,0,0,0,0,
     0,1,0,1,0,0,0,0,
     0,1,0,0,0,0,0,0,
     0,1,0,0,0,0,0,1,
     1,1,1,1,1,1,1,1)

j1<-c(0,0,0,0,1,1,1,1,
    0,0,0,0,0,0,1,0,
    0,0,0,0,0,0,1,0,
    0,0,0,0,0,0,1,0,
    0,0,0,0,0,0,1,0,
    0,0,0,0,0,0,1,0,
    0,0,1,0,0,0,1,0,
    0,0,1,0,0,0,1,0,
    0,0,0,1,1,1,0,0)

k1<-c(1,1,1,0,0,0,1,1,
     0,1,0,0,0,1,0,0,
     0,1,0,0,1,0,0,0,
     0,1,0,1,0,0,0,0,
     0,1,1,0,0,0,0,0,
     0,1,0,1,0,0,0,0,
     0,1,0,0,1,0,0,0,
     0,1,0,0,0,1,0,0,
     1,1,1,0,0,0,1,1)

a2<-c(0,0,0,1,0,0,0,0,
      0,0,0,1,0,0,0,0,
      0,0,0,1,0,0,0,0,
      0,0,1,0,1,0,0,0,
      0,0,1,0,1,0,0,0,
      0,0,1,0,1,0,0,0,
      0,1,1,1,1,1,0,0,
      0,1,0,0,0,1,0,0,
      1,1,1,0,1,1,1,0)

b2<-c(1,1,1,1,1,1,1,0,
      1,0,0,0,0,0,0,1,
      1,0,0,0,0,0,0,1,
      1,0,0,0,0,0,0,1,
      1,1,1,1,1,1,1,0,
      1,0,0,0,0,0,0,1,
      1,0,0,0,0,0,0,1,
      1,0,0,0,0,0,0,1,
      1,1,1,1,1,1,1,0)

c2<-c(0,0,1,1,1,1,0,0,
      0,1,0,0,0,0,1,0,
      1,0,0,0,0,0,0,1,
      1,0,0,0,0,0,0,0,
      1,0,0,0,0,0,0,0,
      1,0,0,0,0,0,0,0,
      1,0,0,0,0,0,0,1,
      0,1,0,0,0,0,1,0,
      0,0,1,1,1,1,0,0)

d2<-c(1,1,1,1,1,1,1,0,
      1,0,0,0,0,0,1,0,
      1,0,0,0,0,0,0,1,
      1,0,0,0,0,0,0,1,
      1,0,0,0,0,0,0,1,
      1,0,0,0,0,0,0,1,
      1,0,0,0,0,0,0,1,
      1,0,0,0,0,0,1,0,
      1,1,1,1,1,1,0,0)

e2<-c(1,1,1,1,1,1,1,1,
      1,0,0,0,0,0,0,1,
      1,0,0,0,0,0,0,0,
      1,0,0,0,0,0,0,0,
      1,1,1,1,1,0,0,0,
      1,0,0,0,0,0,0,0,
      1,0,0,0,0,0,0,0,
      1,0,0,0,0,0,0,0,
      1,1,1,1,1,1,1,1)

j2<-c(0,0,0,0,0,0,1,0,
      0,0,0,0,0,0,1,0,
      0,0,0,0,0,0,1,0,
      0,0,0,0,0,0,1,0,
      0,0,0,0,0,0,1,0,
      0,0,0,0,0,0,1,0,
      0,0,1,0,0,0,1,0,
      0,0,1,0,0,0,1,0,
      0,0,0,1,1,1,0,0)

k2<-c(1,0,0,0,0,1,0,0,
      1,0,0,0,1,0,0,0,
      1,0,0,1,0,0,0,0,
      1,0,1,0,0,0,0,0,
      1,1,0,0,0,0,0,0,
      1,0,1,0,0,0,0,0,
      1,0,0,1,0,0,0,0,
      1,0,0,0,1,0,0,0,
      1,0,0,0,0,1,0,0)

a3<-c(0,0,0,1,0,0,0,0,
      0,0,0,1,0,0,0,0,
      0,0,0,1,0,0,0,0,
      0,0,1,0,1,0,0,0,
      0,0,1,0,1,0,0,0,
      0,1,0,0,0,1,0,0,
      0,1,1,1,1,1,0,0,
      1,0,0,0,0,0,1,0,
      1,1,0,0,0,1,1,0)

b3<-c(1,1,1,1,1,1,1,0,
      0,1,0,0,0,0,0,1,
      0,1,0,0,0,0,0,1,
      0,1,1,1,1,1,1,1,
      0,1,0,0,0,0,0,0,
      0,1,0,0,0,0,0,1,
      0,1,0,0,0,0,0,1,
      0,1,0,0,0,0,0,1,
      1,1,1,1,1,1,1,0)

c3<-c(0,0,1,1,1,1,0,1,
      0,1,0,0,0,0,1,1,
      1,0,0,0,0,0,0,1,
      1,0,0,0,0,0,0,0,
      1,0,0,0,0,0,0,0,
      1,0,0,0,0,0,0,0,
      1,0,0,0,0,0,0,1,
      0,1,0,0,0,0,1,1,
      0,0,1,1,1,1,0,1)

d3<-c(1,1,1,1,1,1,1,0,
      0,1,0,0,0,0,1,0,
      0,1,0,0,0,0,0,1,
      0,1,0,0,0,0,0,1,
      0,1,0,0,0,0,0,1,
      0,1,0,0,0,0,0,1,
      0,1,0,0,0,0,0,1,
      0,1,0,0,0,0,1,0,
      1,1,1,1,1,1,0,0)

e3<-c(1,1,1,1,1,1,1,1,
      0,1,0,0,0,0,0,1,
      0,1,0,0,0,0,0,0,
      0,1,0,0,1,0,0,0,
      0,1,1,1,1,0,0,0,
      0,1,0,0,1,0,0,0,
      0,1,0,0,0,0,0,0,
      0,1,0,0,0,0,0,0,
      1,1,1,1,1,1,1,1)

j3<-c(0,0,0,0,0,1,1,1,
      0,0,0,0,0,0,1,0,
      0,0,0,0,0,0,1,0,
      0,0,0,0,0,0,1,0,
      0,0,0,0,0,0,1,0,
      0,0,0,0,0,0,1,0,
      0,0,1,0,0,0,1,0,
      0,0,1,0,0,0,1,0,
      0,0,0,1,1,1,0,0)

k3<-c(1,1,1,0,0,0,1,0,
      0,1,0,0,0,1,0,0,
      0,1,0,0,1,0,0,0,
      0,1,0,1,0,0,0,0,
      0,1,1,0,0,0,0,0,
      0,1,0,1,0,0,0,0,
      0,1,0,0,1,0,0,0,
      0,1,0,0,0,1,0,0,
      1,1,1,0,0,0,1,0)

aBobot<-c(0,0,0,1,0,0,0,0,
             0,0,0,1,0,0,0,0,
             0,0,0,1,0,0,0,0,
             0,0,1,0,1,0,0,0,
             0,0,1,0,1,0,0,0,
             0,1,0,0,0,1,0,0,
             0,1,1,1,1,1,0,0,
             1,0,0,0,0,0,1,0,
             1,1,0,0,0,1,1,0)

bBobot<-c(1,1,1,1,1,1,1,0,
             0,1,0,0,0,0,0,1,
             0,1,0,0,0,0,0,1,
             0,1,1,1,1,1,1,1,
             0,1,0,0,0,0,0,0,
             0,1,0,0,0,0,0,1,
             0,1,0,0,0,0,0,1,
             0,1,0,0,0,0,0,1,
             1,1,1,1,1,1,1,0)

cBobot<-c(0,0,1,1,1,1,0,1,
             0,1,0,0,0,0,1,1,
             1,0,0,0,0,0,0,1,
             1,0,0,0,0,0,0,0,
             1,0,0,0,0,0,0,0,
             1,0,0,0,0,0,0,0,
             1,0,0,0,0,0,0,1,
             0,1,0,0,0,0,1,1,
             0,0,1,1,1,1,0,1)

dBobot<-c(1,1,1,1,1,1,1,0,
             0,1,0,0,0,0,1,0,
             0,1,0,0,0,0,0,1,
             0,1,0,0,0,0,0,1,
             0,1,0,0,0,0,0,1,
             0,1,0,0,0,0,0,1,
             0,1,0,0,0,0,0,1,
             0,1,0,0,0,0,1,0,
             1,1,1,1,1,1,0,0)

eBobot<-c(1,1,1,1,1,1,1,1,
             0,1,0,0,0,0,0,1,
             0,1,0,0,0,0,0,0,
             0,1,0,0,1,0,0,0,
             0,1,1,1,1,0,0,0,
             0,1,0,0,1,0,0,0,
             0,1,0,0,0,0,0,0,
             0,1,0,0,0,0,0,0,
             1,1,1,1,1,1,1,1)

jBobot<-c(0,0,0,0,0,1,1,1,
             0,0,0,0,0,0,1,0,
             0,0,0,0,0,0,1,0,
             0,0,0,0,0,0,1,0,
             0,0,0,0,0,0,1,0,
             0,0,0,0,0,0,1,0,
             0,0,1,0,0,0,1,0,
             0,0,1,0,0,0,1,0,
             0,0,0,1,1,1,0,0)

kBobot<-c(1,1,1,0,0,0,1,0,
             0,1,0,0,0,1,0,0,
             0,1,0,0,1,0,0,0,
             0,1,0,1,0,0,0,0,
             0,1,1,0,0,0,0,0,
             0,1,0,1,0,0,0,0,
             0,1,0,0,1,0,0,0,
             0,1,0,0,0,1,0,0,
             1,1,1,0,0,0,1,0)

font1 <-matrix(c(a1,b1,c1,d1,e1,j1,k1),ncol=72,byrow=T)
font2 <-matrix(c(a2,b2,c2,d2,e2,j2,k2),ncol=72,byrow=T)
font3 <-matrix(c(a3,b3,c3,d3,e3,j3,k3),ncol=72,byrow=T)
BobotSom<-matrix(c(aBobot,bBobot,cBobot,dBobot,eBobot,jBobot,kBobot),ncol=72,byrow=T)
bobotLvq<-matrix(c(aBobot,bBobot,cBobot,dBobot,eBobot,jBobot,kBobot),ncol=72,byrow=T)
font = array(c(font1,font2,font3),dim=c(7,72,3))


classFeedFoward<-function(font){
    
    #menentukan variabel respn Y
    Y <- c(0,0,0,0,0,0,0)

    cbind(font,Y)

    #menjadikan sebagai bobot inisialisasi untuk layer 1.
    rand_vector <- runif(ncol(font)*nrow(font))

        rand_matrix <- matrix(
        rand_vector,
        nrow = ncol(font),
        ncol = nrow(font),
        byrow = TRUE

    )

    # Buat list yang menyimpan parameter NN yang akan detraining
    my_nn <- list(
    
        # Variabel Input
        input = font,
        # Bobot untuk layer 1
        weights1 = rand_matrix,
        # Bias untuk layer 1
        bias1 = 1.0,
        # Bobot untuk layer 2
        weights2 = matrix(runif(7),ncol=1),
        # Bias untuk layer 2
        bias2 = 1.0,
        # Nilai Aktual
        y = Y,
        # Variabel untuk menyimpan nilai prediksi
        output = matrix(
            rep(0, times = 7),
            ncol = 1
        )
    
    )

    # Fungsi Aktivasi Sigmoid
    sigmoid <- function(x) {
    1.0 / (1.0 + exp(-x))
    }

    # Penjabaran Formulasi Feedforward
    feedforward <- function(nn) {
        nn$layer1 <- sigmoid((nn$input%*%nn$weights1) + nn$bias1)
        nn$output <- sigmoid((nn$layer1%*%nn$weights2) + nn$bias2)
        nn
    }

    my_nn<-feedforward(my_nn)

    data.frame(
        "Predicted" = round(my_nn$output, 3),
        "Actual" = Y
    )

}

classBackProp <- function(my_nn, n){

    loss_function<-function(nn){ 
        sum((nn$y - nn$output) ^ 2) 
    } 

    loss_function(my_nn)

    # Derivasi fungsi aktifasi sigmoid (derivative of sigmoid) 
    sigmoid_derivative <- function(x) { 
        sigmoid(x) * (1.0 - sigmoid(x)) 
    }

    ### Backpropagasi Balik berdasarkan aturan rantai (derivasi) 
    backprop <- function(nn) { 
    
    # d_weights2 adalah aturan rantai (derivasi) bobot 2 
    d_weights2 <- (
        t(nn$layer1) %*% 
        ((nn$y - nn$output) * sigmoid_derivative(nn$output))) 
    
    # d_weights1 adalah aturan rantai (derivasi) bobot 1 dan bobot 2 
    d_weights1 <- ( 
        t(nn$input) %*% 
        (((nn$y - nn$output) * sigmoid_derivative(nn$output))%*% 
            t(nn$weights2)* sigmoid_derivative(nn$layer1))) 
    
        # update bobot menggunakan derivatif(slope) loss function 
        nn$weights1 <- nn$weights1 + d_weights1 
        nn$weights2 <- nn$weights2 + d_weights2 
        nn 
    }

    # Buat data frame untuk menyimpan hasil loss function tiap iterasi 
    loss_df <- data.frame( 
    iteration = 1:n, 
    loss = vector("numeric", length = n))

    # Lakukan pelatihan sebanyak n iterasi kemudian simpan nilai loss 
    for (i in 1:n) { 

        my_nn <- feedforward(my_nn) 
        my_nn <- backprop(my_nn) 
        # Simpan hasil loss function, untuk lihat plot perubahan erorr 
        loss_df$loss[i] <- loss_function(my_nn) 
    
    }


    # lihat hasil prediksi dengan backprop 
    data.frame( 
        "Predicted" = round(my_nn$output, 3), 
        "Actual" = Y 
    )
}

classSOM<-function(font,BobotSom){
    
    data_som<-font
    row.names(data_som)<-c("A1","B1","C1","D1","E1","J1","K1")

    #define BOBOT
    row.names(BobotSom)<-c("A","B","C","D","E","J","K")
    
    n <-nrow(data_som)
    nepoch <- 8
    alpha <- as.vector(rep(0,nepoch))
    final_cluster <-as.vector(rep(0,n))

    for(t in 1:nepoch){
    
    if(t>=1 && t<=4){
        alpha[t]=0.5
    }else{
        alpha[t]=0.5*0.4
    }
    
    for(i in 1:n){

        d1=sum((data_som[i,]-BobotSom["A",])^2)
        d2=sum((data_som[i,]-BobotSom["B",])^2)
        d3=sum((data_som[i,]-BobotSom["C",])^2)
        d4=sum((data_som[i,]-BobotSom["D",])^2)
        d5=sum((data_som[i,]-BobotSom["E",])^2)
        d6=sum((data_som[i,]-BobotSom["J",])^2)
        d7=sum((data_som[i,]-BobotSom["K",])^2)
        
        dsum=data.frame(d1,d2,d3,d4,d5,d6,d7)
        dmin<-min(dsum)
        
        if(BobotSom[i,] == BobotSom[match(dmin,dsum),]){
            BobotSom[i,]<-BobotSom[i,]+alpha[t]*(data_som[i,]-BobotSom[match(dmin,dsum)])
            final_cluster[i] = match(dmin,dsum)
        }else if(BobotSom[i,] != BobotSom[match(dmin,dsum),]){
            BobotSom[i,]<-BobotSom[i,]+alpha[t]*(data_som[i,]-BobotSom[i,])
        }

    }

    # lihat hasil prediksi dengan backprop 
    data.frame( 
        "Predicted" = round(final_cluster),
    )
  }
}

classLvq<-function(font,bobotLvq){

    datasetLvq<-font
    row.names(datasetLvq)<-c("A1","B1","C1","D1","E1","J1","K1")

    #define BOBOT
    row.names(bobotLvq)<-c("A","B","C","D","E","J","K")

    target_lvq<-data.frame(Target=c("A","B","C","D","E","J","K"))

    nLvq <-nrow(datasetLvq)
    nepochLvq <- 10
    alphaLvq <- as.vector(rep(0,nepochLvq))
    final_KelasLvq <-as.vector(rep(0,nLvq))

    for(t in 1:nepochLvq){

    if(t>=1 && t<=4){
        alphaLvq[t]=0.5
    }else{
        alphaLvq[t]=0.5*0.4
    }

    for(i in 1:nLvq){

        d1=sum((datasetLvq[i,]-bobotLvq["A",])^2)
        d2=sum((datasetLvq[i,]-bobotLvq["B",])^2)
        d3=sum((datasetLvq[i,]-bobotLvq["C",])^2)
        d4=sum((datasetLvq[i,]-bobotLvq["D",])^2)
        d5=sum((datasetLvq[i,]-bobotLvq["E",])^2)
        d6=sum((datasetLvq[i,]-bobotLvq["J",])^2)
        d7=sum((datasetLvq[i,]-bobotLvq["K",])^2)

        dsumLvq=data.frame(d1,d2,d3,d4,d5,d6,d7)
        dminLvq<-min(dsumLvq)

            if(bobotLvq[i,] == bobotLvq[match(dminLvq,dsumLvq),]){
                
                if(target_lvq[i,] == "A"){
                bobotLvq[i,]<-bobotLvq[i,]+alpha[t]*(datasetLvq[i,]-bobotLvq[match(dmin,dsum)])
                final_KelasLvq[i] = match(dminLvq,dsumLvq)
                }else{
                bobotLvq["B",]<-bobotLvq["B",]-alphaLvq[t]*(datasetLvq[i,]-bobotLvq["B",])
                bobotLvq["C",]<-bobotLvq["C",]-alphaLvq[t]*(datasetLvq[i,]-bobotLvq["C",])
                bobotLvq["D",]<-bobotLvq["D",]-alphaLvq[t]*(datasetLvq[i,]-bobotLvq["D",])
                bobotLvq["E",]<-bobotLvq["E",]-alphaLvq[t]*(datasetLvq[i,]-bobotLvq["E",])
                bobotLvq["J",]<-bobotLvq["J",]-alphaLvq[t]*(datasetLvq[i,]-bobotLvq["J",])
                bobotLvq["K",]<-bobotLvq["K",]-alphaLvq[t]*(datasetLvq[i,]-bobotLvq["K",])
                }
                
                if(target_lvq[i,] == "B"){
                bobotLvq[i,]<-bobotLvq[i,]+alpha[t]*(datasetLvq[i,]-bobotLvq[match(dmin,dsum)])
                final_KelasLvq[i] = match(dminLvq,dsumLvq)
                }else{
                bobotLvq["A",]<-bobotLvq["A",]-alphaLvq[t]*(datasetLvq[i,]-bobotLvq["A",])
                bobotLvq["C",]<-bobotLvq["C",]-alphaLvq[t]*(datasetLvq[i,]-bobotLvq["C",])
                bobotLvq["D",]<-bobotLvq["D",]-alphaLvq[t]*(datasetLvq[i,]-bobotLvq["D",])
                bobotLvq["E",]<-bobotLvq["E",]-alphaLvq[t]*(datasetLvq[i,]-bobotLvq["E",])
                bobotLvq["J",]<-bobotLvq["J",]-alphaLvq[t]*(datasetLvq[i,]-bobotLvq["J",])
                bobotLvq["K",]<-bobotLvq["K",]-alphaLvq[t]*(datasetLvq[i,]-bobotLvq["K",])
                }

                if(target_lvq[i,] == "C"){
                bobotLvq[i,]<-bobotLvq[i,]+alpha[t]*(datasetLvq[i,]-bobotLvq[match(dmin,dsum)])
                final_KelasLvq[i] = match(dminLvq,dsumLvq)
                }else{
                bobotLvq["A",]<-bobotLvq["A",]-alphaLvq[t]*(datasetLvq[i,]-bobotLvq["A",])
                bobotLvq["B",]<-bobotLvq["B",]-alphaLvq[t]*(datasetLvq[i,]-bobotLvq["B",])
                bobotLvq["D",]<-bobotLvq["D",]-alphaLvq[t]*(datasetLvq[i,]-bobotLvq["D",])
                bobotLvq["E",]<-bobotLvq["E",]-alphaLvq[t]*(datasetLvq[i,]-bobotLvq["E",])
                bobotLvq["J",]<-bobotLvq["J",]-alphaLvq[t]*(datasetLvq[i,]-bobotLvq["J",])
                bobotLvq["K",]<-bobotLvq["K",]-alphaLvq[t]*(datasetLvq[i,]-bobotLvq["K",])
                }
                
                if(target_lvq[i,] == "D"){
                bobotLvq[i,]<-bobotLvq[i,]+alpha[t]*(datasetLvq[i,]-bobotLvq[match(dmin,dsum)])
                final_KelasLvq[i] = match(dminLvq,dsumLvq)
                }else{
                bobotLvq["A",]<-bobotLvq["A",]-alphaLvq[t]*(datasetLvq[i,]-bobotLvq["A",])
                bobotLvq["B",]<-bobotLvq["B",]-alphaLvq[t]*(datasetLvq[i,]-bobotLvq["B",])
                bobotLvq["C",]<-bobotLvq["C",]-alphaLvq[t]*(datasetLvq[i,]-bobotLvq["C",])
                bobotLvq["E",]<-bobotLvq["E",]-alphaLvq[t]*(datasetLvq[i,]-bobotLvq["E",])
                bobotLvq["J",]<-bobotLvq["J",]-alphaLvq[t]*(datasetLvq[i,]-bobotLvq["J",])
                bobotLvq["K",]<-bobotLvq["K",]-alphaLvq[t]*(datasetLvq[i,]-bobotLvq["K",])
                }
                
                if(target_lvq[i,] == "E"){
                bobotLvq[i,]<-bobotLvq[i,]+alpha[t]*(datasetLvq[i,]-bobotLvq[match(dmin,dsum)])
                final_KelasLvq[i] = match(dminLvq,dsumLvq)
                }else{
                bobotLvq["A",]<-bobotLvq["A",]-alphaLvq[t]*(datasetLvq[i,]-bobotLvq["A",])
                bobotLvq["B",]<-bobotLvq["B",]-alphaLvq[t]*(datasetLvq[i,]-bobotLvq["B",])
                bobotLvq["C",]<-bobotLvq["C",]-alphaLvq[t]*(datasetLvq[i,]-bobotLvq["C",])
                bobotLvq["D",]<-bobotLvq["D",]-alphaLvq[t]*(datasetLvq[i,]-bobotLvq["D",])
                bobotLvq["J",]<-bobotLvq["J",]-alphaLvq[t]*(datasetLvq[i,]-bobotLvq["J",])
                bobotLvq["K",]<-bobotLvq["K",]-alphaLvq[t]*(datasetLvq[i,]-bobotLvq["K",])
                }
                
                if(target_lvq[i,] == "J"){
                bobotLvq[i,]<-bobotLvq[i,]+alpha[t]*(datasetLvq[i,]-bobotLvq[match(dmin,dsum)])
                final_KelasLvq[i] = match(dminLvq,dsumLvq)
                }else{
                bobotLvq["A",]<-bobotLvq["A",]-alphaLvq[t]*(datasetLvq[i,]-bobotLvq["A",])
                bobotLvq["B",]<-bobotLvq["B",]-alphaLvq[t]*(datasetLvq[i,]-bobotLvq["B",])
                bobotLvq["C",]<-bobotLvq["C",]-alphaLvq[t]*(datasetLvq[i,]-bobotLvq["C",])
                bobotLvq["D",]<-bobotLvq["D",]-alphaLvq[t]*(datasetLvq[i,]-bobotLvq["D",])
                bobotLvq["E",]<-bobotLvq["E",]-alphaLvq[t]*(datasetLvq[i,]-bobotLvq["E",])
                bobotLvq["K",]<-bobotLvq["K",]-alphaLvq[t]*(datasetLvq[i,]-bobotLvq["K",])
                }
                
                if(target_lvq[i,] == "K"){
                bobotLvq[i,]<-bobotLvq[i,]+alpha[t]*(datasetLvq[i,]-bobotLvq[match(dmin,dsum)])
                final_KelasLvq[i] = match(dminLvq,dsumLvq)
                }else{
                bobotLvq["A",]<-bobotLvq["A",]-alphaLvq[t]*(datasetLvq[i,]-bobotLvq["A",])
                bobotLvq["B",]<-bobotLvq["B",]-alphaLvq[t]*(datasetLvq[i,]-bobotLvq["B",])
                bobotLvq["C",]<-bobotLvq["C",]-alphaLvq[t]*(datasetLvq[i,]-bobotLvq["C",])
                bobotLvq["D",]<-bobotLvq["D",]-alphaLvq[t]*(datasetLvq[i,]-bobotLvq["D",])
                bobotLvq["E",]<-bobotLvq["E",]-alphaLvq[t]*(datasetLvq[i,]-bobotLvq["E",])
                bobotLvq["J",]<-bobotLvq["J",]-alphaLvq[t]*(datasetLvq[i,]-bobotLvq["J",])
                }
            }
        }
    }

    # lihat hasil prediksi dengan backprop 
    data.frame( 
        "Predicted" = round(final_KelasLvq),
    )
}

classHebb<-function(font,training){
    
    datasetHebb<-font
    row.names(datasetHebb)<-c("A1","B1","C1","D1","E1","J1","K1")

    #menjadikan sebagai bobot inisialisasi untuk layer 1.
    rand_vector = runif(ncol(datasetHebb) * datasetHebb)

    rand_matrix = matrix(rand_vector,nrow=7, byrow = T)

    #Bobot untuk layer
    weights1 = matrix(runif(rand_matrix * ncol(datasetHebb)),nrow=7, byrow = T)
              
    #yin
    yin = datasetHebb %*% t(weights1)

    if(training == T){

         #delta
         delta = (t(input)%*%output)

        #weights1
        weights1 = weights1 + t(delta)

    }

    #output
    output = ifelse(yin>0,1,-1)

    # lihat hasil prediksi dengan backprop 
    data.frame( 
        "Predicted" = round(output),
    )

}

classMain<-function(font){
  
  for(i in 1:font){
    
      classFeedFoward(font[,,i])
      classBackProp(classFeedFoward(font[,,i]))
      classSOM(font[,,i],BobotSom)
      classLvq(font[,,i],bobotLvq)
      classHebb(font[,,i],TRUE)
  
  }
  
}
