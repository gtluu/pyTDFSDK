import os
import sqlite3
import pandas as pd
from pyTDFSDK.tims import *
from pyTDFSDK.tsf import *
from pyTDFSDK.util import *
from pyTDFSDK.ctypes_data_structures import *
from pyTDFSDK.error import *


class TsfData(object):
    """
    Class containing metadata from TSF files and methods from TDF-SDK library to work with TSF format data.

    :param bruker_d_folder_name: Path to a Bruker .d directory containing analysis.tsf.
    :type bruker_d_folder_name: str
    :param tdf_sdk: Library initialized by pyTDFSDK.init_tdf_sdk.init_tdf_sdk_api().
    :type tdf_sdk: ctypes.CDLL
    :param use_recalibrated_state: Whether to use recalibrated data (True) or not (False), defaults to True.
    :type use_recalibrated_state: bool
    """
    def __init__(self, bruker_d_folder_name: str, tdf_sdk, use_recalibrated_state=True):
        """
        Constructor Method
        """
        self.api = tdf_sdk
        self.source_file = bruker_d_folder_name
        self.handle = tsf_open(self.api, self.source_file, use_recalibrated_state)
        if self.handle == 0:
            throw_last_tsfdata_error(self.api)
        self.conn = sqlite3.connect(os.path.join(bruker_d_folder_name, 'analysis.tsf'))

        self.analysis = None

        self.get_db_tables()
        self.close_sql_connection()

    def __del__(self):
        """
        Close connection to raw data handle.
        """
        if hasattr(self, 'handle'):
            tsf_close(self.api, self.handle, self.conn)

    def get_db_tables(self):
        """
         Get a dictionary of all tables found in the analysis.tsf SQLite database in which the table names act as keys
         and the tables as a pandas.DataFrame of values; this is stored in pyTDFSDK.classes.TsfData.analysis.
        """
        cursor = self.conn.cursor()
        cursor.execute("SELECT name FROM sqlite_master WHERE type='table';")
        table_names = cursor.fetchall()
        table_names = [table[0] for table in table_names]
        self.analysis = {name: pd.read_sql_query("SELECT * FROM " + name, self.conn) for name in table_names}
        self.analysis['GlobalMetadata'] = {row['Key']: row['Value']
                                           for index, row, in self.analysis['GlobalMetadata'].iterrows()}
        cursor.close()

    def close_sql_connection(self):
        """
        Close the connection to analysis.tsf.
        """
        self.conn.close()


class TdfData(object):
    """
    Class containing metadata from TDF files and methods from TDF-SDK library to work with TDF format data.

    :param bruker_d_folder_name: Path to a Bruker .d directory containing analysis.tdf.
    :type bruker_d_folder_name: str
    :param tdf_sdk: Library initialized by pyTDFSDK.init_tdf_sdk.init_tdf_sdk_api().
    :type tdf_sdk: ctypes.CDLL
    :param use_recalibrated_state: Whether to use recalibrated data (True) or not (False), defaults to True.
    :type use_recalibrated_state: bool
    :param pressure_compensation_strategy: Pressure compensation when opening TDF data
        (pyTDFSDK.ctypes_data_structures.PressureCompensationStrategy.NoPressureCompensation = None,
        pyTDFSDK.ctypes_data_structures.PressureCompensationStrategy.AnalyisGlobalPressureCompensation = analysis
        global pressure compensation,
        pyTDFSDK.ctypes_data_structures.PressureCompensationStrategy.PerFramePressureCompensation = per frame pressure
        compensation), defaults to No Pressure Compensation.
    :type pressure_compensation_strategy: enum.Enum
    """
    def __init__(self, bruker_d_folder_name: str, tdf_sdk, use_recalibrated_state=True,
                 pressure_compensation_strategy=PressureCompensationStrategy.NoPressureCompensation):
        """
        Constructor Method
        """
        self.api = tdf_sdk
        self.source_file = bruker_d_folder_name
        self.handle = tims_open_v2(self.api, self.source_file, pressure_compensation_strategy, use_recalibrated_state)
        if self.handle == 0:
            throw_last_timsdata_error(self.api)
        self.conn = sqlite3.connect(os.path.join(bruker_d_folder_name, 'analysis.tdf'))

        self.analysis = None

        self.get_db_tables()
        self.close_sql_connection()

    def __del__(self):
        """
        Close connection to raw data handle.
        """
        if hasattr(self, 'handle'):
            tims_close(self.api, self.handle, self.conn)

    def get_db_tables(self):
        """
         Get a dictionary of all tables found in the analysis.tsf SQLite database in which the table names act as keys
         and the tables as a pandas.DataFrame of values; this is stored in pyTDFSDK.classes.TsfData.analysis.
        """
        cursor = self.conn.cursor()
        cursor.execute("SELECT name FROM sqlite_master WHERE type='table';")
        table_names = cursor.fetchall()
        table_names = [table[0] for table in table_names]
        self.analysis = {name: pd.read_sql_query("SELECT * FROM " + name, self.conn) for name in table_names}
        self.analysis['GlobalMetadata'] = {row['Key']: row['Value']
                                           for index, row, in self.analysis['GlobalMetadata'].iterrows()}
        cursor.close()

    def close_sql_connection(self):
        """
        Close the connection to analysis.tdf.
        """
        self.conn.close()
