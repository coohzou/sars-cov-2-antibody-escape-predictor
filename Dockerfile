FROM python:3.12-slim

WORKDIR /app

RUN apt-get update && apt-get install -y --no-install-recommends \
    gcc \
    && rm -rf /var/lib/apt/lists/*

COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt gunicorn

COPY app.py .
COPY templates/ templates/
COPY utils/ utils/
COPY data/ data/

ENV PORT=7860
EXPOSE 7860

CMD gunicorn --bind 0.0.0.0:${PORT} --timeout 300 --workers 1 --threads 1 app:app
