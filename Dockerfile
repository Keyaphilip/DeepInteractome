# DeepInteractome Dockerfile
FROM python:3.10-slim

# Set environment variables
ENV PYTHONDONTWRITEBYTECODE 1
ENV PYTHONUNBUFFERED 1
ENV PYTHONPATH /app

# Set work directory
WORKDIR /app

# Install system dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    && rm -rf /var/lib/apt/lists/*

# Install python dependencies (filtering out heavy packages like torch for lighter builds if needed, 
# but models require torch and joblib so we install from requirements.txt)
COPY requirements.txt /app/
RUN pip install --upgrade pip && \
    pip install -r requirements.txt --no-cache-dir

# Copy project
COPY . /app/

# Expose port
EXPOSE 8000

# Run the FastAPI application
CMD ["uvicorn", "api.main:app", "--host", "0.0.0.0", "--port", "8000"]
