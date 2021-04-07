try:
    import loguru
except:
    print("Package loguru not found. Install loguru to get log file.")

logger = loguru.logger

logger.add("CactusTool.log", format="{time:YYYY-MM-DD at HH:mm:ss} | {level} | {message}", level="INFO", rotation="10 MB", backtrace=True, diagnose=True)