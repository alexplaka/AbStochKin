#  Copyright (c) 2025, Alex Plakantonakis.
#
#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.

import pathlib
import logging
from datetime import datetime


class CustomFormatter(logging.Formatter):
    """ Format the specified log record by processing its name attribute. """

    def format(self, record: logging.LogRecord) -> str:
        """
        This method removes the 'logging_config.' substring and the '.py' suffix
        from the record's name attribute before formatting it using the parent
        class's format method.

        Parameters
        ----------
        record : logging.LogRecord
            The log record to be formatted.

        Returns
        -------
        str
            The formatted log message.
        """
        record.name = record.name.replace('logging_config.', '').replace('.py', '')
        return super().format(record)


logger = logging.getLogger(__name__)  # Create a logger

logger.setLevel(logging.DEBUG)  # Set the logging level

# Create a file handler and a stream handler
pathlib.Path("logs").mkdir(parents=True, exist_ok=True)
file_handler = logging.FileHandler(f"logs/{datetime.now().strftime("%Y-%m-%d_%H:%M:%S")}.log")
file_handler.setLevel(logging.DEBUG)

stream_handler = logging.StreamHandler()
stream_handler.setLevel(logging.WARNING)

# Create a formatter and add it to the handlers
formatter = CustomFormatter("{asctime} - {name} - {levelname} - {message}",
                            datefmt='%Y-%m-%d %H:%M:%S',
                            style="{")

file_handler.setFormatter(formatter)
stream_handler.setFormatter(formatter)

# Add the handlers to the logger
logger.addHandler(file_handler)
logger.addHandler(stream_handler)

logger.info('Log initiated.')  # Log first message
