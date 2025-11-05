python -m venv venv
.\venv\Scripts\activate

pip install numpy==2.3.4
pip install torch==2.4.1 torchvision==0.19.1 torchaudio==2.4.1 --index-url https://download.pytorch.org/whl/cu124
pip install transformers==4.47.1
pip install requests==2.32.5
pip install chromadb==0.5.20

pip install hf_xet