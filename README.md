# DeepT-Hybrid: Deep Learning and Hybrid Features for Protein Classification

This project presents a machine learning framework to classify protein sequences as "Positive" or "Negative". The project utilizes a hybrid feature extraction method that combines two different approaches:

1.  **Deep Features:** Semantic and structural features (embeddings) are extracted from protein sequences using **ProtBERT-BFD**, a powerful protein language model.
2.  **Chemical Features:** 365 different physicochemical features are calculated for the proteins, including amino acid composition, hydrophobicity, and polarity, based on the iTTCA-RF method.

These two feature sets are combined ("hybrid features") to train and evaluate several classification models, including Convolutional Neural Networks (CNN), XGBoost, and Random Forest.

---

## 🚀 Models and Results

Multiple models were trained and tested to find the best-performing approach. The final results on the test set are as follows:

| Model Name | Features Used | Test Accuracy |
| :------------------------------ | :--------------------------------- | :-------------: |
| iTTCA-RF (Random Forest) | Chemical Features Only | ~86% |
| DeepT-i (Inception CNN) | ProtBERT Deep Features Only | ~87% |
| DeepT-Hybrid-Inception (CNN) | **Hybrid** (ProtBERT + Chemical) | ~88% |
| **DeepT-Hybrid-XGB (XGBoost)** | **Hybrid** (ProtBERT + Chemical) | **~90.7%** |

The best performance was achieved by the **DeepT-Hybrid-XGB** model, which leverages the combined power of deep and chemical features with an XGBoost classifier.

---

## Installation & Requirements

Before running the project, ensure you have the following libraries installed. You can create a `requirements.txt` file and install them using the command `pip install -r requirements.txt`.

**requirements.txt**
```
tensorflow
torch
transformers
xgboost
scikit-learn
pandas
numpy
matplotlib
shap
keras-tuner
visualkeras
pydot
graphviz
```

Alternatively, you can install the libraries individually:
```bash
pip install tensorflow torch transformers xgboost scikit-learn pandas numpy matplotlib shap keras-tuner visualkeras pydot graphviz
```

---

## Required Files

To train and test the models, the following **3 files** must be present in the project's root directory (the same location as the `.ipynb` files):

1.  **`TrainSet.txt`**: A text file in a FASTA-like format containing the protein sequences for training. Each sequence header must include a `Positive` or `Negative` label.
2.  **`TestSet.txt`**: A text file containing the protein sequences for testing.
3.  **`PAAC.txt`**: A file containing the "Pseudo Amino Acid Composition" (PAAC) parameters used during the chemical feature extraction process.

---

## Usage

The project is provided in three Jupyter Notebook files:

1.  **`DeepT_Hybrid.ipynb` (Recommended)**:
    * This is the main notebook containing all the experiments mentioned above. It allows you to run and compare every model from a single file.

2.  **`DeepT_Hybrid_Inception.ipynb`**:
    * This notebook focuses specifically on the Inception-based CNN models (`DeepT-i` and `DeepT-Hybrid-Inception`).

3.  **`DeepT_Hybrid_XGB.ipynb`**:
    * This notebook focuses on the XGBoost model (`DeepT-Hybrid-XGB`) and includes SHAP analysis for feature importance.

**Steps:**

1.  Install the required libraries listed above.
2.  Upload `TrainSet.txt`, `TestSet.txt` files, which stay in the Human Dataset 2 folder and `PAAC.txt` file in the Google Colab environment. 
3.  Open your preferred notebook (we recommend starting with `DeepT_Hybrid.ipynb`) in a Jupyter or Google Colab environment.
4.  Execute all cells sequentially (`Run All`).
5.  The results (accuracy, AUC score, classification reports) will be displayed in the output of the corresponding cells.

**Example File Structure:**
```
.
├── DeepT_Hybrid.ipynb
├── DeepT_Hybrid_Inception.ipynb
├── DeepT_Hybrid_XGB.ipynb
├── TrainSet.txt
├── TestSet.txt
├── PAAC.txt
└── README.md
```

---

## 💻 Technologies Used

* Python 3
* TensorFlow & Keras (for building and training deep learning models)
* Keras-Tuner (for hyperparameter optimization)
* Hugging Face Transformers (for integrating pre-trained models like ProtBERT)
* XGBoost & Scikit-learn (for traditional ML and model evaluation metrics)
* Pandas & NumPy (for data processing and numerical computations)
* Matplotlib (for basic plotting and data visualization)
* visualkeras (for visualizing Keras model architectures)
* pydot & graphviz (backend support for visualization tools)
* SHAP (for explainability of model predictions)
