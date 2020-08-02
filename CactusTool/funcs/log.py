import loguru

logger = loguru.logger

logger.add("CactusTool.log", format="{time:YYYY-MM-DD at HH:mm:ss} | {level} | {message}", level="INFO", rotation="10 MB", backtrace=True, diagnose=True)