FROM python:3.12-slim

WORKDIR /app

COPY pyproject.toml README.md ./
COPY extract ./extract
COPY validation ./validation
COPY orchestration ./orchestration
COPY config ./config
COPY dbt ./dbt

RUN python -m pip install --upgrade pip && pip install -e ".[extract,dbt]"

# Pipeline entrypoint is wired in across later milestones.
CMD ["python", "-c", "print('SARVault scaffold ready')"]
