# TLSformer

TLSformer is an R package that can assist immune researchers in identifying which cells/spots belong to TLS regions or in determining the relative distance between TLS and non-TLS regions for cells/spots. It achieves this through using spatial transcriptomics data to predict single-cell RNA-seq data or other spatial transcriptomics data. TLSformer has advantages in that it does not require a large amount of training data and can perform zero-shot TLS prediction for new data from cancer types it was not directly trained on.

If you want to deep into the TLSformer you can find the manuscript in XXX and we also deposit a tutorial in the XXX.
- [Manuscript](https://www.example.com)
- [Tutorial](https://www.example.com)

## Getting Started

This is a mini example shows how to use TLSformer to predict TLS.

### Installing

Before you use TLSformer, you can firstly create a conda environment.

Using yaml file to create TLSformer conda environment

    conda env create -f tlsformer_env.yml

And after successfully creating the environment, you can find the python of this environment like in this path

    /home/xfan/miniconda3/envs/TLSformer_env/bin/python

This path will be the finally used python environment path, and then download the 10x Visium breast cancer pre-trained gene word encoder and demo data in the Google Cloud. The saved path of this pre-train gene word encoder will be used in the next work flow.
- [pre-trained gene word encoder](https://drive.google.com/drive/folders/1qLsl22T3IU2EEyXYM3z52_8MLNsFDyjO?usp=drive_link)
- [demo data](https://drive.google.com/drive/folders/1DZJ-f_RjpnRUszXNKm_KRGXpbHcwsEBK?usp=drive_link)

Install TLSformer by devtools in R

    devtools::install_github("Jinglab/TLSformer")
    
Alternatively, you can download the [TLSformer_1.0.tar.gz](https://github.com/Jinglab/TLSformer/blob/main/TLSformer_1.0.tar.gz) file from this GitHub repository and install it locally.

    install.packages("~/MLTLS_package/TLSformer_1.0.tar.gz")
    
### Run TLSformer 

1.Load package and demo data

    library(Seurat)
    library(TLSformer)
    library(reticulate)
    library(tidyverse)
    st_dat_train <- readRDS("~/MLTLS_package/demo_data/bc_st_demo_data.RDS")
    st_dat_pred <- readRDS("~/MLTLS_package/demo_data/melanoma_st_demo_data.rds")

2.Set parameters

    sen_len = 260
    save_inseu = TRUE
    genes_representor = "~/MLTLS_package/demo_data/pretrained_models_rank260/genelist.txt"
    envir_path = "/home/xfan/miniconda3/envs/tlsformer_env/bin/python"
    pretrained_model = "TLSformer_BERT"
    pretrained_model_path = "~/MLTLS_package/demo_data/pretrained_models_rank260/"
    save_checkpoint_path = "~/MLTLS_package/demo_data/"
    batch_size = 2
    train_K = 5
    train_Q = 5
    train_episodes = 300

3.Generate sentences
    
    # Training data
    st_dat_train <- generate_sentences(
      seu_obj = st_dat_train,
      sen_len = sen_len,
      region_info = st_dat_train@meta.data$region,
      save_inseu = save_inseu,
      genes_representor = genes_representor,
      envir_path = envir_path
    )
    
    # Predicton data
    st_dat_pred <- generate_sentences(
      seu_obj = st_dat_pred,
      sen_len = sen_len,
      region_info = st_dat_pred@meta.data$region,
      save_inseu = save_inseu,
      genes_representor = genes_representor,
      envir_path = envir_path
    )

4.Training TLSformer
    
    # Training
    st_dat_train <- run_tlsformer_train(
        seu_obj = st_dat_train,
        pretrained_model = pretrained_model,
        sen_len = sen_len,
        pretrained_model_path = pretrained_model_path,
        save_checkpoint_path = save_checkpoint_path,
        batch_size = batch_size,
        train_K = train_K,
        train_Q = train_Q,
        train_episodes = train_episodes,
        envir_path = envir_path
    )

5.Use trained TLSformer to predict

    # Run prediction
    st_dat_pred <- run_tlsformer_pred(
                        seu_obj = st_dat_pred,
                        pretrained_model_path = pretrained_model_path,
                        save_checkpoint_path = save_checkpoint_path,
                        envir_path = envir_path,
                        pretrained_model = pretrained_model,
                        sen_len=sen_len)
    # Normalization -- 0-1 scale
    st_dat_pred$relative_distance <- 1- (st_dat_pred$relative_distance - min(st_dat_pred$relative_distance))/(max(st_dat_pred$relative_distance)-min(st_dat_pred$relative_distance))
    SpatialFeaturePlot(st_dat_pred,features = c("region","relative_distance"))

## Built With
  - [Python](https://www.python.org/) 
  - [R](https://www.contributor-covenant.org/](https://www.r-project.org/about.html)) 

## Lab website

  - **YingJing Lab** - *Guangzhou, China* - (https://www.yingjinglab.com)
