# TLSformer v1.0.0 </a>

### Xiaokai Fan<sup></sup>,  Mei Xie<sup></sup>, Xinying Xue*, Xiaohui Fan*

TLSformer is a computational tool that can assist immune researchers for identifying which TLS cells/spots or calculate determining the probability of a cell/spot belong to TLS regions. Different from previous methods, TLSformer leverages bidirectional encoder representation from transformer (BERT) and meta-learning to acquire general knowledge from the limited available TLS-related information. This knowledge is then transferred to identify single cells or spots within TLS regions from scRNA-seq data or spatial transcriptomics data. 

## Requirements and Installation
This toolkit is written in both R and Python programming languages. The core BERT and meta-learning algorithm are implemented in Python, while the initial data preparation and functions usage are written in R.

### Installation of TLSformer

To use TLSformer, firstly create a conda environment.

Using [yaml file](https://github.com/Jinglab/TLSformer/blob/main/tlsformer_env.yml) to create TLSformer conda environment

    conda env create -f tlsformer_env.yml

After successfully creating the environment, the python path of this environment can be found like this

    /home/xfan/miniconda3/envs/TLSformer_env/bin/python

This path will be the finally used python environment path, and then download the 10x Visium breast cancer pre-trained BERT and demo data in the Google cloud. The saved location of this pre-trained BERT will be utilized in the following workflow steps.
- [breast cancer pre-trained gene word encoder](https://drive.google.com/drive/folders/1qLsl22T3IU2EEyXYM3z52_8MLNsFDyjO?usp=drive_link)
- [demo data](https://drive.google.com/drive/folders/1DZJ-f_RjpnRUszXNKm_KRGXpbHcwsEBK?usp=drive_link)

Install TLSformer by devtools in R
[![R >4.0](https://img.shields.io/badge/R-%3E%3D4.0-brightgreen)](https://www.r-project.org/)

    devtools::install_github("Jinglab/TLSformer")
    
Alternatively, you can download the [TLSformer_1.0.tar.gz](https://github.com/Jinglab/TLSformer/blob/main/TLSformer_1.0.tar.gz) file from this GitHub repository and install it locally.

    install.packages("home/xfan/MLTLS_package/TLSformer_1.0.tar.gz")

## Quick Start

### Run TLSformer 
To use TLSformer, we require formatted `.csv` files as input (i.e. read in by pandas). 

1.Load package and demo data

    library(Seurat)
    library(TLSformer)
    library(reticulate)
    library(tidyverse)
    st_dat_train <- readRDS("~/MLTLS_package/demo_data/bc_st_demo_data.rds")
    st_dat_pred <- readRDS("~/MLTLS_package/demo_data/melanoma_st_demo_data.rds")

2.Set parameters

    sen_len = 260
    save_inseu = TRUE
    genes_representor = "~/MLTLS_package/demo_data/pretrained_models_rank260/genelist.txt"
    envir_path = "/home/xfan/miniconda3/envs/TLSformer_env/bin/python"
    pretrained_model = "TLSformer_BERT"
    pretrained_model_path = "~/MLTLS_package/demo_data/pretrained_models_rank260/"
    save_checkpoint_path = "~/MLTLS_package/demo_data/"
    batch_size = 1
    train_K = 2
    train_Q = 2
    train_episodes = 600

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

    # Prediction
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
  - [PyTorch](https://pytorch.org/)
  - [R](https://www.contributor-covenant.org/](https://www.r-project.org/about.html))


## Lab website

  - **YingJing Lab** - *Guangzhou, China* - https://www.yingjinglab.com
