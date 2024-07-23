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
        cursor.execute("SELECT name FROM sqlite_master WHERE type='view';")
        view_names = cursor.fetchall()
        view_names = [view[0] for view in view_names]
        names = table_names + view_names
        self.analysis = {name: pd.read_sql_query("SELECT * FROM " + name, self.conn) for name in names}
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
        compensation), defaults to Global Pressure Compensation.
    :type pressure_compensation_strategy: enum.Enum
    """
    def __init__(self, bruker_d_folder_name: str, tdf_sdk, use_recalibrated_state=True,
                 pressure_compensation_strategy=PressureCompensationStrategy.AnalyisGlobalPressureCompensation):
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
        cursor.execute("SELECT name FROM sqlite_master WHERE type='view';")
        view_names = cursor.fetchall()
        view_names = [view[0] for view in view_names]
        names = table_names + view_names
        self.analysis = {name: pd.read_sql_query("SELECT * FROM " + name, self.conn) for name in names}
        self.analysis['GlobalMetadata'] = {row['Key']: row['Value']
                                           for index, row, in self.analysis['GlobalMetadata'].iterrows()}
        cursor.close()

    def close_sql_connection(self):
        """
        Close the connection to analysis.tdf.
        """
        self.conn.close()


class TsfSpectrum(object):
    """
    Class for parsing and storing spectrum metadata and data arrays from BAF format data.

    :param tsf_data: TsfData object containing metadata from analysis.tsf database.
    :type tsf_data: pyTDFSDK.classes.TsfData
    :param frame: ID of the frame of interest.
    :type frame: int
    :param mode: Data array mode, either "profile", "centroid", or "raw".
    :type mode: str
    :param profile_bins: Number of bins to bin spectrum to.
    :type profile_bins: int
    :param mz_encoding: m/z encoding command line parameter, either "64" or "32".
    :type mz_encoding: int
    :param intensity_encoding: Intensity encoding command line parameter, either "64" or "32".
    :type intensity_encoding: int
    """

    def __init__(self, tsf_data, frame: int, mode: str, profile_bins=0, mz_encoding=64, intensity_encoding=64):
        """
        Constructor Method
        """
        self.tsf_data = tsf_data
        self.scan_number = None
        self.scan_type = None
        self.ms_level = None
        self.mz_array = None
        self.intensity_array = None
        self.mobility_array = None
        self.polarity = None
        self.centroided = None
        self.retention_time = None
        self.coord = None
        self.total_ion_current = None
        self.base_peak_mz = None
        self.base_peak_intensity = None
        self.high_mz = None
        self.low_mz = None
        self.target_mz = None
        self.isolation_lower_offset = None
        self.isolation_upper_offset = None
        self.selected_ion_mz = None
        self.selected_ion_intensity = None
        self.selected_ion_mobility = None
        self.selected_ion_ccs = None
        self.charge_state = None
        self.activation = None
        self.collision_energy = None
        self.frame = frame
        self.parent_frame = None
        self.parent_scan = None
        self.ms2_no_precursor = False
        self.mode = mode
        self.profile_bins = profile_bins
        self.mz_encoding = mz_encoding
        self.intensity_encoding = intensity_encoding

        if 'MaldiApplicationType' not in self.tsf_data.analysis['GlobalMetadata'].keys():
            self.get_lcms_tsf_data()
        elif 'MaldiApplicationType' in self.tsf_data.analysis['GlobalMetadata'].keys():
            self.get_maldi_tsf_data()

    def get_lcms_tsf_data(self):
        frames_dict = self.tsf_data.analysis['Frames'][self.tsf_data.analysis['Frames']['Id'] ==
                                                       self.frame].to_dict(orient='records')[0]
        self.polarity = frames_dict['Polarity']
        self.centroided = get_centroid_status(self.mode)[0]
        self.retention_time = float(frames_dict['Time']) / 60
        self.mz_array, self.intensity_array = extract_tsf_spectrum(self.tsf_data,
                                                                   self.frame,
                                                                   self.mode,
                                                                   self.profile_bins,
                                                                   self.mz_encoding,
                                                                   self.intensity_encoding)
        if self.mz_array is not None and self.intensity_array is not None and \
                self.mz_array.size != 0 and self.intensity_array.size != 0 and \
                self.mz_array.size == self.intensity_array.size:
            self.total_ion_current = sum(self.intensity_array)
            base_peak_index = np.where(self.intensity_array == np.max(self.intensity_array))
            self.base_peak_mz = self.mz_array[base_peak_index][0].astype(float)
            self.base_peak_intensity = self.intensity_array[base_peak_index][0].astype(float)
            self.high_mz = float(max(self.mz_array))
            self.low_mz = float(min(self.mz_array))
            # MS1
            if int(frames_dict['MsMsType']) == 0:
                self.scan_type = 'MS1 spectrum'
                self.ms_level = 1
            elif int(frames_dict['MsMsType']) in [2, 8, 9]:
                framemsmsinfo_dict = self.tsf_data.analysis['FrameMsMsInfo'][self.tsf_data.analysis['FrameMsMsInfo']['Frame'] ==
                                                                             self.frame].to_dict(orient='records')[0]
                # Auto MS/MS
                if int(frames_dict['ScanMode']) == 1:
                    self.scan_type = 'MSn spectrum'
                    self.ms_level = 2
                    self.target_mz = float(framemsmsinfo_dict['TriggerMass'])
                    self.isolation_lower_offset = float(framemsmsinfo_dict['IsolationWidth']) / 2
                    self.isolation_upper_offset = float(framemsmsinfo_dict['IsolationWidth']) / 2
                    self.selected_ion_mz = float(framemsmsinfo_dict['TriggerMass'])
                    self.charge_state = framemsmsinfo_dict['PrecursorCharge']
                    self.activation = 'collision-induced dissociation'
                    self.collision_energy = float(framemsmsinfo_dict['CollisionEnergy'])
                    self.parent_frame = int(framemsmsinfo_dict['Parent'])
                # MRM
                elif int(frames_dict['ScanMode']) == 2:
                    self.scan_type = 'MSn spectrum'
                    self.ms_level = 2
                    self.target_mz = float(framemsmsinfo_dict['TriggerMass'])
                    self.isolation_lower_offset = float(framemsmsinfo_dict['IsolationWidth']) / 2
                    self.isolation_upper_offset = float(framemsmsinfo_dict['IsolationWidth']) / 2
                    self.selected_ion_mz = float(framemsmsinfo_dict['TriggerMass'])
                    self.charge_state = framemsmsinfo_dict['PrecursorCharge']
                    self.activation = 'collision-induced dissociation'
                    self.collision_energy = float(framemsmsinfo_dict['CollisionEnergy'])
                # bbCID
                elif int(frames_dict['ScanMode']) == 4:
                    self.scan_type = 'MSn spectrum'
                    self.ms_level = 2
                    self.activation = 'collision-induced dissociation'
                    self.collision_energy = float(framemsmsinfo_dict['CollisionEnergy'])
                    self.ms2_no_precursor = True

    def get_maldi_tsf_data(self):
        frames_dict = self.tsf_data.analysis['Frames'][self.tsf_data.analysis['Frames']['Id'] ==
                                                       self.frame].to_dict(orient='records')[0]
        maldiframeinfo_dict = self.tsf_data.analysis['MaldiFrameInfo'][self.tsf_data.analysis['MaldiFrameInfo']['Frame'] ==
                                                                       self.frame].to_dict(orient='records')[0]
        self.coord = get_maldi_coords(self.tsf_data, maldiframeinfo_dict)
        self.polarity = frames_dict['Polarity']
        self.centroided = get_centroid_status(self.mode)[0]
        self.retention_time = 0
        self.mz_array, self.intensity_array = extract_tsf_spectrum(self.tsf_data,
                                                                   self.frame,
                                                                   self.mode,
                                                                   self.profile_bins,
                                                                   self.mz_encoding,
                                                                   self.intensity_encoding)
        if self.mz_array is not None and self.intensity_array is not None and \
                self.mz_array.size != 0 and self.intensity_array.size != 0 and \
                self.mz_array.size == self.intensity_array.size:
            self.total_ion_current = sum(self.intensity_array)
            base_peak_index = np.where(self.intensity_array == np.max(self.intensity_array))
            self.base_peak_mz = self.mz_array[base_peak_index][0].astype(float)
            self.base_peak_intensity = self.intensity_array[base_peak_index][0].astype(float)
            self.high_mz = float(max(self.mz_array))
            self.low_mz = float(min(self.mz_array))
            # MS1
            if int(frames_dict['MsMsType']) == 0:
                self.scan_type = 'MS1 spectrum'
                self.ms_level = 1
            # MS/MS
            elif int(frames_dict['MsMsType']) in [2, 8, 9]:
                msms_mode_id = self.tsf_data.analysis['PropertyDefinitions'][(self.tsf_data.analysis['PropertyDefinitions']['PermanentName'] ==
                                                                             'Mode_ScanMode')].to_dict(orient='records')[0]['Id']
                msms_mode = self.tsf_data.analysis['Properties'][(self.tsf_data.analysis['Properties']['Frame'] == self.frame) &
                                                                 (self.tsf_data.analysis['Properties']['Property'] == msms_mode_id)].to_dict(orient='records')[0]['Value']
                framemsmsinfo_dict = self.tsf_data.analysis['FrameMsMsInfo'][self.tsf_data.analysis['FrameMsMsInfo']['Frame'] ==
                                                                             self.frame].to_dict(orient='records')[0]
                # MALDI MS/MS, coded as MRM in the schema
                if msms_mode == 3:
                    self.scan_type = 'MSn spectrum'
                    self.ms_level = 2
                    self.target_mz = float(framemsmsinfo_dict['TriggerMass'])
                    self.isolation_lower_offset = float(framemsmsinfo_dict['IsolationWidth']) / 2
                    self.isolation_upper_offset = float(framemsmsinfo_dict['IsolationWidth']) / 2
                    self.selected_ion_mz = float(framemsmsinfo_dict['TriggerMass'])
                    self.charge_state = framemsmsinfo_dict['PrecursorCharge']
                    self.activation = 'collision-induced dissociation'
                    self.collision_energy = float(framemsmsinfo_dict['CollisionEnergy'])
                # MALDI bbCID
                elif msms_mode == 5:
                    self.scan_type = 'MSn spectrum'
                    self.ms_level = 2
                    self.activation = 'collision-induced dissociation'
                    self.collision_energy = float(framemsmsinfo_dict['CollisionEnergy'])
                    self.ms2_no_precursor = True


class TdfSpectrum(object):
    """
    Class for parsing and storing spectrum metadata and data arrays from BAF format data.

    :param tdf_data: TdfData object containing metadata from analysis.tdf database.
    :type tdf_data: pyTDFSDK.classes.TdfData
    :param frame: ID of the frame of interest.
    :type frame: int
    :param mode: Data array mode, either "profile", "centroid", or "raw".
    :type mode: str
    :param precursor: ID of the precursor of interest for ddaPASEF data. If specified, overrides the frame ID during
        data parsing.
    :type precursor: int
    :param diapasef_window: Dictionary containing a row of metadata from the
        pyTDFSDK.classesTdfData.analysis['DiaFrameMsMsWindows'] table required for parsing diaPASEF data.
    :type diapasef_window: dict
    :param profile_bins: Number of bins to bin spectrum to.
    :type profile_bins: int
    :param mz_encoding: m/z encoding command line parameter, either "64" or "32".
    :type mz_encoding: int
    :param intensity_encoding: Intensity encoding command line parameter, either "64" or "32".
    :type intensity_encoding: int
    :param mobility_encoding: Mobility encoding command line parameter, either "64" or "32".
    :type mobility_encoding: int
    """

    def __init__(self, tdf_data, frame: int, mode: str, precursor=0, diapasef_window=None, profile_bins=0,
                 mz_encoding=64, intensity_encoding=64, mobility_encoding=64, exclude_mobility=False):
        """
        Constructor Method
        """
        self.tdf_data = tdf_data
        self.scan_number = None
        self.scan_type = None
        self.ms_level = None
        self.mz_array = None
        self.intensity_array = None
        self.mobility_array = None
        self.polarity = None
        self.centroided = None
        self.retention_time = None
        self.coord = None
        self.total_ion_current = None
        self.base_peak_mz = None
        self.base_peak_intensity = None
        self.high_mz = None
        self.low_mz = None
        self.target_mz = None
        self.isolation_lower_offset = None
        self.isolation_upper_offset = None
        self.selected_ion_mz = None
        self.selected_ion_intensity = None
        self.selected_ion_mobility = None
        self.selected_ion_ccs = None
        self.charge_state = None
        self.activation = None
        self.collision_energy = None
        self.frame = frame
        self.parent_frame = None
        self.parent_scan = None
        self.ms2_no_precursor = False
        self.precursor = precursor
        self.diapasef_window = diapasef_window  # Should be a row from TdfData.analysis['DiaFrameMsMsWindows'] table.
        self.mode = mode
        self.profile_bins = profile_bins
        self.mz_encoding = mz_encoding
        self.intensity_encoding = intensity_encoding
        self.mobility_encoding = mobility_encoding
        self.exclude_mobility = exclude_mobility

        if self.precursor != 0:
            self.frame = 0

        if 'MaldiApplicationType' not in self.tdf_data.analysis['GlobalMetadata'].keys():
            if self.frame != 0 and self.precursor == 0:
                self.get_lcms_tdf_data()
            elif self.frame == 0 and self.precursor != 0:
                self.get_ddapasef_precursor_data()
        elif 'MaldiApplicationType' in self.tdf_data.analysis['GlobalMetadata'].keys():
            self.get_maldi_tdf_data()

    def get_lcms_tdf_data(self):
        frames_dict = self.tdf_data.analysis['Frames'][self.tdf_data.analysis['Frames']['Id'] ==
                                                       self.frame].to_dict(orient='records')[0]
        self.polarity = frames_dict['Polarity']
        self.centroided = get_centroid_status(self.mode, self.exclude_mobility)[0]
        self.retention_time = float(frames_dict['Time']) / 60
        # MS1
        if int(frames_dict['MsMsType']) == 0:
            self.scan_type = 'MS1 spectrum'
            self.ms_level = 1
            if not self.exclude_mobility:
                self.mz_array, self.intensity_array, self.mobility_array = extract_3d_tdf_spectrum(self.tdf_data,
                                                                                                   self.frame,
                                                                                                   0,
                                                                                                   int(frames_dict['NumScans']))
            elif self.exclude_mobility:
                self.mz_array, self.intensity_array = extract_2d_tdf_spectrum(self.tdf_data,
                                                                              self.frame,
                                                                              0,
                                                                              int(frames_dict['NumScans']),
                                                                              self.mode,
                                                                              self.profile_bins,
                                                                              self.mz_encoding,
                                                                              self.intensity_encoding)
            if self.mz_array is not None and self.intensity_array is not None and \
                    self.mz_array.size != 0 and self.intensity_array.size != 0 and \
                    self.mz_array.size == self.intensity_array.size:
                self.total_ion_current = sum(self.intensity_array)
                base_peak_index = np.where(self.intensity_array == np.max(self.intensity_array))
                self.base_peak_mz = self.mz_array[base_peak_index][0].astype(float)
                self.base_peak_intensity = self.intensity_array[base_peak_index][0].astype(float)
                self.high_mz = float(max(self.mz_array))
                self.low_mz = float(min(self.mz_array))
        # diaPASEF
        elif int(frames_dict['ScanMode']) == 9 and int(frames_dict['MsMsType']) == 9:
            if self.diapasef_window is not None:
                if not self.exclude_mobility:
                    self.mz_array, self.intensity_array, self.mobility_array = extract_3d_tdf_spectrum(self.tdf_data,
                                                                                                       self.frame,
                                                                                                       int(self.diapasef_window['ScanNumBegin']),
                                                                                                       int(self.diapasef_window['ScanNumEnd']))
                elif self.exclude_mobility:
                    self.mz_array, self.intensity_array = extract_2d_tdf_spectrum(self.tdf_data,
                                                                                  self.frame,
                                                                                  int(self.diapasef_window['ScanNumBegin']),
                                                                                  int(self.diapasef_window['ScanNumEnd']),
                                                                                  self.mode,
                                                                                  self.profile_bins,
                                                                                  self.mz_encoding,
                                                                                  self.intensity_encoding)
                if self.mz_array is not None and self.intensity_array is not None and \
                        self.mz_array.size != 0 and self.intensity_array.size != 0 and \
                        self.mz_array.size == self.intensity_array.size:
                    self.scan_type = 'MSn spectrum'
                    self.ms_level = 2
                    self.target_mz = float(self.diapasef_window['IsolationMz'])
                    self.isolation_lower_offset = float(self.diapasef_window['IsolationWidth']) / 2
                    self.isolation_upper_offset = float(self.diapasef_window['IsolationWidth']) / 2
                    self.selected_ion_mz = float(self.diapasef_window['IsolationMz'])
                    self.collision_energy = self.diapasef_window['CollisionEnergy']
                    self.total_ion_current = sum(self.intensity_array)
                    base_peak_index = np.where(self.intensity_array == np.max(self.intensity_array))
                    self.base_peak_mz = self.mz_array[base_peak_index][0].astype(float)
                    self.base_peak_intensity = self.intensity_array[base_peak_index][0].astype(float)
                    self.high_mz = float(max(self.mz_array))
                    self.low_mz = float(min(self.mz_array))
            else:
                raise Exception("diapasef_window from TdfData.analysis['DiaFrameMsMsWindows'] table required to parse "
                                "diaPASEF data.")
        # bbCID
        elif int(frames_dict['ScanMode']) == 4 and int(frames_dict['MsMsType']) == 2:
            framemsmsinfo_dict = self.tdf_data.analysis['FrameMsMsInfo'][self.tdf_data.analysis['FrameMsMsInfo']['Frame'] ==
                                                                         self.frame].to_dict(orient='records')[0]
            if not self.exclude_mobility:
                self.mz_array, self.intensity_array, self.mobility_array = extract_3d_tdf_spectrum(self.tdf_data,
                                                                                                   self.frame,
                                                                                                   0,
                                                                                                   int(frames_dict['NumScans']))
            elif self.exclude_mobility:
                self.mz_array, self.intensity_array = extract_2d_tdf_spectrum(self.tdf_data,
                                                                              self.frame,
                                                                              0,
                                                                              int(frames_dict['NumScans']),
                                                                              self.mode,
                                                                              self.profile_bins,
                                                                              self.mz_encoding,
                                                                              self.intensity_encoding)
            if self.mz_array is not None and self.intensity_array is not None and \
                    self.mz_array.size != 0 and self.intensity_array.size != 0 and \
                    self.mz_array.size == self.intensity_array.size:
                self.scan_type = 'MSn spectrum'
                self.ms_level = 2
                self.activation = 'collision-induced dissociation'
                self.collision_energy = float(framemsmsinfo_dict['CollisionEnergy'])
                self.ms2_no_precursor = True
                self.total_ion_current = sum(self.intensity_array)
                base_peak_index = np.where(self.intensity_array == np.max(self.intensity_array))
                self.base_peak_mz = self.mz_array[base_peak_index][0].astype(float)
                self.base_peak_intensity = self.intensity_array[base_peak_index][0].astype(float)
                self.high_mz = float(max(self.mz_array))
                self.low_mz = float(min(self.mz_array))
        # MRM
        elif int(frames_dict['ScanMode']) == 2 and int(frames_dict['MsMsType']) == 2:
            framemsmsinfo_dict = self.tdf_data.analysis['FrameMsMsInfo'][self.tdf_data.analysis['FrameMsMsInfo']['Frame'] ==
                                                                         self.frame].to_dict(orient='records')[0]
            self.mz_array, self.intensity_array = extract_2d_tdf_spectrum(self.tdf_data,
                                                                          self.frame,
                                                                          0,
                                                                          int(frames_dict['NumScans']),
                                                                          self.mode,
                                                                          self.profile_bins,
                                                                          self.mz_encoding,
                                                                          self.intensity_encoding)
            if self.mz_array is not None and self.intensity_array is not None and \
                    self.mz_array.size != 0 and self.intensity_array.size != 0 and \
                    self.mz_array.size == self.intensity_array.size:
                self.scan_type = 'MSn spectrum'
                self.ms_level = 2
                self.target_mz = float(framemsmsinfo_dict['TriggerMass'])
                self.isolation_lower_offset = float(framemsmsinfo_dict['IsolationWidth']) / 2
                self.isolation_upper_offset = float(framemsmsinfo_dict['IsolationWidth']) / 2
                self.selected_ion_mz = float(framemsmsinfo_dict['TriggerMass'])
                self.charge_state = framemsmsinfo_dict['PrecursorCharge']
                self.activation = 'collision-induced dissociation'
                self.collision_energy = framemsmsinfo_dict['CollisionEnergy']
                self.total_ion_current = sum(self.intensity_array)
                base_peak_index = np.where(self.intensity_array == np.max(self.intensity_array))
                self.base_peak_mz = self.mz_array[base_peak_index][0].astype(float)
                self.base_peak_intensity = self.intensity_array[base_peak_index][0].astype(float)
                self.high_mz = float(max(self.mz_array))
                self.low_mz = float(min(self.mz_array))
        # prmPASEF
        elif int(frames_dict['ScanMode']) == 10 and int(frames_dict['MsMsType']) == 10:
            prmframemsmsinfo_dict = self.tdf_data.analysis['PrmFrameMsMsInfo'][self.tdf_data.analysis['PrmFrameMsMsInfo']['Frame'] ==
                                                                               self.frame].to_dict(orient='records')[0]
            prmtargets_dict = self.tdf_data.analysis['PrmTargets'][self.tdf_data.analysis['PrmTargets']['Id'] ==
                                                                   int(prmframemsmsinfo_dict['Target'])].to_dict(orient='records')[0]
            self.mz_array, self.intensity_array = extract_2d_tdf_spectrum(self.tdf_data,
                                                                          self.frame,
                                                                          int(prmframemsmsinfo_dict['ScanNumBegin']),
                                                                          int(prmframemsmsinfo_dict['ScanNumEnd']),
                                                                          self.mode,
                                                                          self.profile_bins,
                                                                          self.mz_encoding,
                                                                          self.intensity_encoding)
            if self.mz_array is not None and self.intensity_array is not None and \
                    self.mz_array.size != 0 and self.intensity_array.size != 0 and \
                    self.mz_array.size == self.intensity_array.size:
                self.scan_type = 'MSn spectrum'
                self.ms_level = 2
                self.target_mz = float(prmframemsmsinfo_dict['IsolationMz'])
                self.isolation_lower_offset = float(prmframemsmsinfo_dict['IsolationWidth']) / 2
                self.isolation_upper_offset = float(prmframemsmsinfo_dict['IsolationWidth']) / 2
                self.selected_ion_mz = float(prmframemsmsinfo_dict['IsolationMz'])
                self.selected_ion_mobility = float(prmtargets_dict['OneOverK0'])
                self.charge_state = prmtargets_dict['Charge']
                self.activation = 'collision-induced dissociation'
                self.collision_energy = prmframemsmsinfo_dict['CollisionEnergy']
                if not np.isnan(prmtargets_dict['Charge']):
                    self.selected_ion_ccs = tims_oneoverk0_to_ccs_for_mz(self.tdf_data.api,
                                                                         float(prmtargets_dict['OneOverK0']),
                                                                         int(prmtargets_dict['Charge']),
                                                                         float(prmframemsmsinfo_dict['IsolationMz']))
                self.total_ion_current = sum(self.intensity_array)
                base_peak_index = np.where(self.intensity_array == np.max(self.intensity_array))
                self.base_peak_mz = self.mz_array[base_peak_index][0].astype(float)
                self.base_peak_intensity = self.intensity_array[base_peak_index][0].astype(float)
                self.high_mz = float(max(self.mz_array))
                self.low_mz = float(min(self.mz_array))

    def get_ddapasef_precursor_data(self):
        precursor_dict = self.tdf_data.analysis['Precursors'][self.tdf_data.analysis['Precursors']['Id'] ==
                                                              self.precursor].to_dict(orient='records')[0]
        frames_dict = self.tdf_data.analysis['Frames'][self.tdf_data.analysis['Frames']['Id'] ==
                                                       precursor_dict['Parent']].to_dict(orient='records')[0]
        pasefframemsmsinfo_dicts = self.tdf_data.analysis['PasefFrameMsMsInfo'][self.tdf_data.analysis['PasefFrameMsMsInfo']['Precursor'] ==
                                                                                self.precursor].to_dict(orient='records')
        self.polarity = frames_dict['Polarity']
        self.centroided = get_centroid_status(self.mode, self.exclude_mobility)[0]
        self.retention_time = float(frames_dict['Time']) / 60
        self.mz_array, self.intensity_array = extract_ddapasef_precursor_spectrum(self.tdf_data,
                                                                                  pasefframemsmsinfo_dicts,
                                                                                  self.mode,
                                                                                  self.profile_bins,
                                                                                  self.mz_encoding,
                                                                                  self.intensity_encoding)
        if self.mz_array is not None and self.intensity_array is not None and \
                self.mz_array.size != 0 and self.intensity_array.size != 0 and \
                self.mz_array.size == self.intensity_array.size:
            self.scan_type = 'MSn spectrum'
            self.ms_level = 2
            self.target_mz = float(precursor_dict['AverageMz'])
            self.isolation_lower_offset = float(pasefframemsmsinfo_dicts[0]['IsolationWidth']) / 2
            self.isolation_upper_offset = float(pasefframemsmsinfo_dicts[0]['IsolationWidth']) / 2
            self.selected_ion_mz = float(precursor_dict['LargestPeakMz'])
            self.selected_ion_intensity = float(precursor_dict['Intensity'])
            self.selected_ion_mobility = tims_scannum_to_oneoverk0(self.tdf_data.api,
                                                                   self.tdf_data.handle,
                                                                   int(precursor_dict['Parent']),
                                                                   np.array([int(precursor_dict['ScanNumber'])]))[0]
            self.charge_state = precursor_dict['Charge']
            self.activation = 'collision-induced dissociation'
            self.collision_energy = pasefframemsmsinfo_dicts[0]['CollisionEnergy']
            self.parent_frame = int(precursor_dict['Parent'])
            self.parent_scan = int(precursor_dict['ScanNumber'])
            if not np.isnan(precursor_dict['Charge']):
                self.selected_ion_ccs = tims_oneoverk0_to_ccs_for_mz(self.tdf_data.api,
                                                                     self.selected_ion_mobility,
                                                                     int(precursor_dict['Charge']),
                                                                     float(precursor_dict['LargestPeakMz']))
            self.total_ion_current = sum(self.intensity_array)
            base_peak_index = np.where(self.intensity_array == np.max(self.intensity_array))
            self.base_peak_mz = self.mz_array[base_peak_index][0].astype(float)
            self.base_peak_intensity = self.intensity_array[base_peak_index][0].astype(float)
            self.high_mz = float(max(self.mz_array))
            self.low_mz = float(min(self.mz_array))

    def get_maldi_tdf_data(self):
        frames_dict = self.tdf_data.analysis['Frames'][self.tdf_data.analysis['Frames']['Id'] ==
                                                       self.frame].to_dict(orient='records')[0]
        maldiframeinfo_dict = self.tdf_data.analysis['MaldiFrameInfo'][self.tdf_data.analysis['MaldiFrameInfo']['Frame'] ==
                                                                       self.frame].to_dict(orient='records')[0]
        self.coord = get_maldi_coords(self.tdf_data, maldiframeinfo_dict)
        self.polarity = frames_dict['Polarity']
        self.centroided = get_centroid_status(self.mode)[0]
        self.retention_time = 0
        # MS1
        if int(frames_dict['MsMsType']) == 0:
            self.scan_type = 'MS1 spectrum'
            self.ms_level = 1
            if not self.exclude_mobility:
                self.mz_array, self.intensity_array, self.mobility_array = extract_3d_tdf_spectrum(self.tdf_data,
                                                                                                   self.frame,
                                                                                                   0,
                                                                                                   int(frames_dict['NumScans']))
            elif self.exclude_mobility:
                self.mz_array, self.intensity_array = extract_2d_tdf_spectrum(self.tdf_data,
                                                                              self.frame,
                                                                              0,
                                                                              int(frames_dict['NumScans']),
                                                                              self.mode,
                                                                              self.profile_bins,
                                                                              self.mz_encoding,
                                                                              self.intensity_encoding)
            if self.mz_array is not None and self.intensity_array is not None and \
                    self.mz_array.size != 0 and self.intensity_array.size != 0 and \
                    self.mz_array.size == self.intensity_array.size:
                self.total_ion_current = sum(self.intensity_array)
                base_peak_index = np.where(self.intensity_array == np.max(self.intensity_array))
                self.base_peak_mz = self.mz_array[base_peak_index][0].astype(float)
                self.base_peak_intensity = self.intensity_array[base_peak_index][0].astype(float)
                self.high_mz = float(max(self.mz_array))
                self.low_mz = float(min(self.mz_array))
        elif int(frames_dict['MsMsType']) in [2, 8, 9]:
            msms_mode_id = self.tdf_data.analysis['PropertyDefinitions'][self.tdf_data.analysis['PropertyDefinitions']['PermanentName'] ==
                                                                         'Mode_ScanMode'].to_dict(orient='records')[0]['Id']
            msms_mode = self.tdf_data.analysis['Properties'][(self.tdf_data.analysis['Properties']['Frame'] == self.frame) &
                                                             (self.tdf_data.analysis['Properties']['Property'] == msms_mode_id)].to_dict(orient='records')[0]['Value']
            # MALDI MS/MS, coded as MRM in the schema
            if msms_mode == 3:
                framemsmsinfo_dict = self.tdf_data.analysis['FrameMsMsInfo'][self.tdf_data.analysis['FrameMsMsInfo']['Frame'] ==
                                                                             maldiframeinfo_dict['Frame']].to_dict(orient='records')[0]
                if not self.exclude_mobility:
                    self.mz_array, self.intensity_array, self.mobility_array = extract_3d_tdf_spectrum(self.tdf_data,
                                                                                                       self.frame,
                                                                                                       0,
                                                                                                       int(frames_dict['NumScans']))
                elif self.exclude_mobility:
                    self.mz_array, self.intensity_array = extract_2d_tdf_spectrum(self.tdf_data,
                                                                                  self.frame,
                                                                                  0,
                                                                                  int(frames_dict['NumScans']),
                                                                                  self.mode,
                                                                                  self.profile_bins,
                                                                                  self.mz_encoding,
                                                                                  self.intensity_encoding)
                if self.mz_array is not None and self.intensity_array is not None and \
                        self.mz_array.size != 0 and self.intensity_array.size != 0 and \
                        self.mz_array.size == self.intensity_array.size:
                    self.scan_type = 'MSn spectrum'
                    self.ms_level = 2
                    self.target_mz = float(framemsmsinfo_dict['TriggerMass'])
                    self.isolation_lower_offset = float(framemsmsinfo_dict['IsolationWidth']) / 2
                    self.isolation_upper_offset = float(framemsmsinfo_dict['IsolationWidth']) / 2
                    self.selected_ion_mz = float(framemsmsinfo_dict['TriggerMass'])
                    self.charge_state = framemsmsinfo_dict['PrecursorCharge']
                    self.activation = 'collision-induced dissociation'
                    self.collision_energy = float(framemsmsinfo_dict['CollisionEnergy'])
                    self.total_ion_current = sum(self.intensity_array)
                    base_peak_index = np.where(self.intensity_array == np.max(self.intensity_array))
                    self.base_peak_mz = self.mz_array[base_peak_index][0].astype(float)
                    self.base_peak_intensity = self.intensity_array[base_peak_index][0].astype(float)
                    self.high_mz = float(max(self.mz_array))
                    self.low_mz = float(min(self.mz_array))
            # MALDI bbCID
            elif msms_mode == 5:
                framemsmsinfo_dict = self.tdf_data.analysis['FrameMsMsInfo'][self.tdf_data.analysis['FrameMsMsInfo']['Frame'] ==
                                                                             maldiframeinfo_dict['Frame']].to_dict(orient='records')[0]
                if not self.exclude_mobility:
                    self.mz_array, self.intensity_array, self.mobility_array = extract_3d_tdf_spectrum(self.tdf_data,
                                                                                                       self.frame,
                                                                                                       0,
                                                                                                       int(frames_dict['NumScans']))
                elif self.exclude_mobility:
                    self.mz_array, self.intensity_array = extract_2d_tdf_spectrum(self.tdf_data,
                                                                                  self.frame,
                                                                                  0,
                                                                                  int(frames_dict['NumScans']),
                                                                                  self.mode,
                                                                                  self.profile_bins,
                                                                                  self.mz_encoding,
                                                                                  self.intensity_encoding)
                if self.mz_array is not None and self.intensity_array is not None and \
                        self.mz_array.size != 0 and self.intensity_array.size != 0 and \
                        self.mz_array.size == self.intensity_array.size:
                    self.scan_type = 'MSn spectrum'
                    self.ms_level = 2
                    self.activation = 'collision-induced dissociation'
                    self.collision_energy = float(framemsmsinfo_dict['CollisionEnergy'])
                    self.ms2_no_precursor = True
                    self.total_ion_current = sum(self.intensity_array)
                    base_peak_index = np.where(self.intensity_array == np.max(self.intensity_array))
                    self.base_peak_mz = self.mz_array[base_peak_index][0].astype(float)
                    self.base_peak_intensity = self.intensity_array[base_peak_index][0].astype(float)
                    self.high_mz = float(max(self.mz_array))
                    self.low_mz = float(min(self.mz_array))
            # MALDI prmPASEF
            elif msms_mode == 12:
                if self.diapasef_window is not None:
                    if not self.exclude_mobility:
                        self.mz_array, self.intensity_array, self.mobility_array = extract_3d_tdf_spectrum(self.tdf_data,
                                                                                                           self.frame,
                                                                                                           int(self.diapasef_window['ScanNumBegin']),
                                                                                                           int(self.diapasef_window['ScanNumEnd']))
                    elif self.exclude_mobility:
                        self.mz_array, self.intensity_array = extract_2d_tdf_spectrum(self.tdf_data,
                                                                                      self.frame,
                                                                                      int(self.diapasef_window['ScanNumBegin']),
                                                                                      int(self.diapasef_window['ScanNumEnd']),
                                                                                      self.mode,
                                                                                      self.profile_bins,
                                                                                      self.mz_encoding,
                                                                                      self.intensity_encoding)
                    if self.mz_array is not None and self.intensity_array is not None and \
                            self.mz_array.size != 0 and self.intensity_array.size != 0 and \
                            self.mz_array.size == self.intensity_array.size:
                        self.scan_type = 'MSn spectrum'
                        self.ms_level = 2
                        self.target_mz = float(self.diapasef_window['IsolationMz'])
                        self.isolation_lower_offset = float(self.diapasef_window['IsolationWidth']) / 2
                        self.isolation_upper_offset = float(self.diapasef_window['IsolationWidth']) / 2
                        self.selected_ion_mz = float(self.diapasef_window['IsolationMz'])
                        self.collision_energy = self.diapasef_window['CollisionEnergy']
                        self.total_ion_current = sum(self.intensity_array)
                        base_peak_index = np.where(self.intensity_array == np.max(self.intensity_array))
                        self.base_peak_mz = self.mz_array[base_peak_index][0].astype(float)
                        self.base_peak_intensity = self.intensity_array[base_peak_index][0].astype(float)
                        self.high_mz = float(max(self.mz_array))
                        self.low_mz = float(min(self.mz_array))
                else:
                    raise Exception("diapasef_window from TdfData.analysis['DiaFrameMsMsWindows'] table required to "
                                    "parse diaPASEF data.")
