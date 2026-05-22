"""
Compound flood return-period analysis (Tokyo).

Estimates bivariate rainfall–tide return periods with copula-based joint distributions.
See ../README.md for methodology and usage.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
from scipy import stats
import scipy.optimize as optimize
from scipy.stats import norm
from copulae import GaussianCopula, StudentCopula, ClaytonCopula, GumbelCopula, FrankCopula
import warnings
warnings.filterwarnings('ignore')
from pathlib import Path

DATA_DIR = Path(__file__).resolve().parent / "data"


# Set the default font to Arial for all plots
plt.rcParams['font.family'] = 'Arial'

# Marginal candidates: Hong Kong paper set plus log-normal (Tokyo extension).
MARGINAL_DISTRIBUTIONS = (
    'gamma', 'genextreme', 'genpareto', 'logistic', 'norm', 'gumbel_r', 'nakagami', 'invgauss',
    'lognorm',
)

def resolve_marginal_candidates(excluded_distributions=None):
    """Return standard marginal candidates, optionally omitting excluded forms."""
    if not excluded_distributions:
        return MARGINAL_DISTRIBUTIONS
    excluded = set(excluded_distributions)
    return tuple(d for d in MARGINAL_DISTRIBUTIONS if d not in excluded)

class SimpleJoeCopula:
    """Internal helper."""
    def __init__(self, dim=2):
        self.dim = dim
        self.params = 1.5
    
    def fit(self, data):
        """Internal helper."""
        def neg_log_likelihood(theta):
            if theta < 1:
                return 1e10
            
            theta = float(theta)
            log_likes = np.zeros(len(data))
            
            for i in range(len(data)):
                u_i, v_i = data[i]
                pdf_val = self._pdf_calculate(u_i, v_i, theta)
                log_likes[i] = np.log(pdf_val)
            
            return -np.sum(log_likes)
        
        initial_thetas = [1.5, 2.0, 3.0, 5.0]
        
        best_theta = None
        best_nll = float('inf')
        
        for init_t in initial_thetas:
            try:
                result = optimize.minimize(
                    neg_log_likelihood, 
                    x0=init_t, 
                    bounds=[(1.0, None)],
                    method='L-BFGS-B'
                )
                
                if result.success and result.fun < best_nll:
                    best_nll = result.fun
                    best_theta = result.x[0]
            except:
                continue
        
        if best_theta is None:
            best_theta = 1.5
        
        self.params = float(best_theta)
        return self
    
    def cdf(self, u_v):
        """Internal helper."""
        if isinstance(u_v, list):
            u, v = u_v
        else:
            u, v = u_v[0], u_v[1]
        
        return self._cdf_calculate(u, v, self.params)
    
    def pdf(self, u_v):
        """Internal helper."""
        if isinstance(u_v, list) or (isinstance(u_v, np.ndarray) and u_v.ndim == 1):
            u, v = u_v[0], u_v[1]
        else:
            if u_v.ndim == 2:
                return np.array([self.pdf(point) for point in u_v])
            else:
                u, v = u_v[:, 0], u_v[:, 1]
        
        return self._pdf_calculate(u, v, self.params)
    
    def _cdf_calculate(self, u, v, theta):
        """Internal helper."""
        # C_θ(u, v) = 1 - [(1-u)^θ + (1-v)^θ - (1-u)^θ(1-v)^θ]^{1/θ}
        term1 = (1 - u) ** theta
        term2 = (1 - v) ** theta
        return 1 - ((term1 + term2 - term1 * term2) ** (1/theta))
    
    def _pdf_calculate(self, u, v, theta):
        """Internal helper."""
        term_u = (1 - u) ** (theta - 1)
        term_v = (1 - v) ** (theta - 1)
        
        term1 = (1 - u) ** theta
        term2 = (1 - v) ** theta
        
        sum_term = term1 + term2 - term1 * term2
        
        pdf = (theta - 1) * term_u * term_v * (sum_term ** (1/theta - 2)) * (1 + (term1 * term2 / sum_term))
        
        return pdf

# Datum correction functions
def load_datum_correction_file(offset_file_path):
    """Internal helper."""
    try:
        if not offset_file_path or not pd.io.common.file_exists(offset_file_path):
            # logging omitted
            return None
            
        offset_data = pd.read_csv(offset_file_path, encoding='utf-8')
        # logging omitted
        
        required_cols = ['Year', 'Month', 'Day', 'offset']
        missing_cols = [col for col in required_cols if col not in offset_data.columns]
        if missing_cols:
            # logging omitted
            return None
        
        for col in ['Year', 'Month', 'Day']:
            offset_data[col] = pd.to_numeric(offset_data[col], errors='coerce')
        offset_data['offset'] = pd.to_numeric(offset_data['offset'], errors='coerce')
        
        original_rows = len(offset_data)
        offset_data = offset_data.dropna(subset=required_cols)
        if len(offset_data) < original_rows:
            # logging omitted
            pass
        
        offset_data = create_date_column(offset_data, "offset data")
        if offset_data is None or offset_data.empty:
            # logging omitted
            return None
        
        offset_data = offset_data.sort_values('DATE').reset_index(drop=True)
        
        # logging omitted
        # logging omitted
        
        return offset_data
        
    except Exception as e:
        # logging omitted
        return None

def apply_datum_correction(tide_data, offset_data, sea_level_col='Sea_Level'):
    """Internal helper."""
    if offset_data is None or offset_data.empty:
        # logging omitted
        return tide_data.copy()
    
    corrected_tide_data = tide_data.copy()
    
    if not pd.api.types.is_datetime64_any_dtype(corrected_tide_data['DATE']):
        corrected_tide_data['DATE'] = pd.to_datetime(corrected_tide_data['DATE'], errors='coerce')
    if not pd.api.types.is_datetime64_any_dtype(offset_data['DATE']):
        offset_data['DATE'] = pd.to_datetime(offset_data['DATE'], errors='coerce')
    
    offset_data_sorted = offset_data.sort_values('DATE').reset_index(drop=True)
    
    corrected_tide_data['Applied_Offset'] = 0.0
    corrected_tide_data[f'{sea_level_col}_Original'] = corrected_tide_data[sea_level_col].copy()
    corrected_tide_data['Offset_Period'] = ''
    
    correction_count = 0
    
    # logging omitted
    for i, row in offset_data_sorted.iterrows():
        # logging omitted
    
        pass
    earliest_offset_row = offset_data_sorted.iloc[0]
    earliest_offset_date = earliest_offset_row['DATE']
    earliest_offset_value = earliest_offset_row['offset']
    
    # logging omitted
    
    for idx, tide_row in corrected_tide_data.iterrows():
        tide_date = tide_row['DATE']
        
        applicable_offset = None
        period_description = ""
        
        if tide_date < earliest_offset_date:
            applicable_offset = earliest_offset_value
        else:
            for i, offset_row in offset_data_sorted.iterrows():
                current_date = offset_row['DATE']
                
                if i == len(offset_data_sorted) - 1:
                    if tide_date >= current_date:
                        applicable_offset = offset_row['offset']
                        break
                else:
                    next_date = offset_data_sorted.iloc[i+1]['DATE']
                    if current_date <= tide_date < next_date:
                        applicable_offset = offset_row['offset']
                        break
        
        if applicable_offset is not None:
            corrected_tide_data.loc[idx, sea_level_col] = tide_row[sea_level_col] + applicable_offset
            corrected_tide_data.loc[idx, 'Applied_Offset'] = applicable_offset
            corrected_tide_data.loc[idx, 'Offset_Period'] = period_description
            correction_count += 1
        else:
            corrected_tide_data.loc[idx, 'Applied_Offset'] = 0.0
    
    early_data_count = len(corrected_tide_data[corrected_tide_data['DATE'] < earliest_offset_date])
    early_corrected_count = len(corrected_tide_data[
        (corrected_tide_data['DATE'] < earliest_offset_date) & 
        (corrected_tide_data['Applied_Offset'] == earliest_offset_value)
    ])
    
    # logging omitted
    # logging omitted
    # logging omitted
    # logging omitted
    # logging omitted
    # logging omitted
    # logging omitted
    if early_data_count > 0:
        # logging omitted
    
        pass
    if correction_count > 0:
        original_mean = corrected_tide_data[f'{sea_level_col}_Original'].mean()
        corrected_mean = corrected_tide_data[sea_level_col].mean()
        mean_offset = corrected_tide_data['Applied_Offset'].mean()
        
        # logging omitted
        # logging omitted
        # logging omitted
        # logging omitted
        
        offset_usage = corrected_tide_data['Applied_Offset'].value_counts().sort_index()
        # logging omitted
        for offset_val, count in offset_usage.items():
            # logging omitted
        
            pass
        period_usage = corrected_tide_data['Offset_Period'].value_counts()
        # logging omitted
        for period, count in period_usage.items():
                # logging omitted
    
            pass
    return corrected_tide_data

def export_corrected_tide_data(corrected_tide_data, station_name, sea_level_col='Sea_Level'):
    """Internal helper."""
    try:
        export_data = corrected_tide_data.copy()
        
        final_export_data = export_data[['Year', 'Month', 'Day', sea_level_col]].copy()
        
        final_export_data = final_export_data.rename(columns={sea_level_col: 'Sea level(m)'})
        
        final_export_data['Year'] = final_export_data['Year'].astype(int)
        final_export_data['Month'] = final_export_data['Month'].astype(int) 
        final_export_data['Day'] = final_export_data['Day'].astype(int)
        
        final_export_data['Sea level(m)'] = final_export_data['Sea level(m)'].round(4)
        
        output_filename = f'{station_name}-tide-correction.csv'
        
        final_export_data.to_csv(output_filename, index=False, encoding='utf-8-sig')
        
        if 'Applied_Offset' in export_data.columns:
            applied_corrections = (export_data['Applied_Offset'] != 0).sum()
            correction_stats = export_data['Applied_Offset'].describe()
            
            # logging omitted
            # logging omitted
            # logging omitted
            # logging omitted
            # logging omitted
            # logging omitted
            # logging omitted
            # logging omitted
        else:
            # logging omitted
            # logging omitted
            # logging omitted
            # logging omitted
        
            pass
        return output_filename
        
    except Exception as e:
        # logging omitted
        return None

# Load data
def load_data(tide_data_path="WAG_tide.csv", rainfall_data_paths=None, export_rainfall=True, 
              offset_file_path=None):
    """
    Load tide and rainfall data with support for multiple rainfall sources and datum correction.

    Parameters:
    -----------
    tide_data_path : str
        Path to the tide data CSV file
    rainfall_data_paths : list of str or str
        Paths to rainfall data files. Can be a single string or a list.
        Earlier files in the list have higher priority for duplicate dates.
    export_rainfall : bool
        Whether to export the merged rainfall data to CSV
    offset_file_path : str, optional
        Path to the datum correction file containing Year, Month, Day, offset columns.
        If provided, datum correction will be applied to tide data.

    Returns:
    --------
    tuple
        (tide_data, rainfall_data) as pandas DataFrames
        Note: tide_data will have datum correction applied if offset_file_path is provided
    """
    try:
        # Handle input formats
        if rainfall_data_paths is None:
            # Default to a common pattern if none provided
            rainfall_data_paths = ["daily_WGL_RF_ALL.csv"] # Adjust default if needed
            print(f"Warning: No rainfall paths provided, using default: {rainfall_data_paths}")
        elif isinstance(rainfall_data_paths, str):
            rainfall_data_paths = [rainfall_data_paths]

        # Load tide data
        tide_data = pd.read_csv(tide_data_path, encoding='utf-8')
        print(f"Original tide data shape from '{tide_data_path}': {tide_data.shape}")

        # Process tide data date columns
        date_cols_tide = ['Year', 'Month', 'Day']
        missing_date_cols_tide = [col for col in date_cols_tide if col not in tide_data.columns]
        if missing_date_cols_tide:
            print(f"Error: Missing date columns in tide data: {missing_date_cols_tide}")
            return None, None

        # Convert date columns to numeric, coercing errors
        for col in date_cols_tide:
            tide_data[col] = pd.to_numeric(tide_data[col], errors='coerce')

        original_tide_rows_date = len(tide_data)
        tide_data = tide_data.dropna(subset=date_cols_tide)
        rows_dropped_tide = original_tide_rows_date - len(tide_data)
        if rows_dropped_tide > 0:
            print(f"Warning: Removed {rows_dropped_tide} rows from tide data with invalid numeric dates")

        # Convert valid date columns to integer type
        for col in date_cols_tide:
             tide_data[col] = tide_data[col].astype(int)

        # Create DATE column for tide data
        tide_data = create_date_column(tide_data, "tide")
        if tide_data is None:
            print("Error creating DATE column for tide data.")
            return None, None

        # Process sea level column
        sea_level_col_original = 'Sea level(m)' # Original name
        sea_level_col_target = 'Sea_Level'    # Target name

        if sea_level_col_original in tide_data.columns:
            print(f"Renaming tide column '{sea_level_col_original}' to '{sea_level_col_target}'.")
            tide_data.rename(columns={sea_level_col_original: sea_level_col_target}, inplace=True)
        elif sea_level_col_target not in tide_data.columns:
            # Check for common alternative names
            alt_names = ['SeaLevel', 'WaterLevel', 'Water Level(m)', 'Sea_level(m)']
            found_alt = None
            for alt in alt_names:
                if alt in tide_data.columns:
                    print(f"Warning: Target column '{sea_level_col_target}' not found. Found alternative '{alt}' and renaming it.")
                    tide_data.rename(columns={alt: sea_level_col_target}, inplace=True)
                    found_alt = True
                    break
            if not found_alt:
                print(f"Error: Target sea level column '{sea_level_col_target}' and common alternatives not found.")
                return None, None

        # Convert sea level to numeric, coercing errors
        tide_data[sea_level_col_target] = pd.to_numeric(tide_data[sea_level_col_target], errors='coerce')
        missing_sea_level = tide_data[sea_level_col_target].isna().sum()
        if missing_sea_level > 0:
            print(f"Warning: Removing {missing_sea_level} tide rows with missing/non-numeric sea level.")
            tide_data = tide_data.dropna(subset=[sea_level_col_target])

        # Apply datum correction if offset file is provided
        if offset_file_path:
            # logging omitted
            # logging omitted
            offset_data = load_datum_correction_file(offset_file_path)
            if offset_data is not None:
                tide_data = apply_datum_correction(tide_data, offset_data, sea_level_col_target)
                # logging omitted
            else:
                # logging omitted
                pass
        else:
            # logging omitted

            pass
        # Load and combine multiple rainfall files
        combined_rainfall = process_multiple_rainfall_files(rainfall_data_paths)
        if combined_rainfall is None or combined_rainfall.empty:
            print("Error: Failed to process rainfall data.")
            return None, None

        # Export combined rainfall data if requested
        if export_rainfall and combined_rainfall is not None and not combined_rainfall.empty:
            export_merged_rainfall(combined_rainfall, tide_data_path)

        # --- Final Checks ---
        if combined_rainfall.empty or tide_data.empty:
            print("Error: Tide or Rainfall DataFrame is empty after initial loading/cleaning.")
            return None, None

        # Ensure DATE columns are datetime objects
        if not pd.api.types.is_datetime64_any_dtype(tide_data['DATE']):
             tide_data['DATE'] = pd.to_datetime(tide_data['DATE'], errors='coerce')
             tide_data = tide_data.dropna(subset=['DATE'])
        if not pd.api.types.is_datetime64_any_dtype(combined_rainfall['DATE']):
             combined_rainfall['DATE'] = pd.to_datetime(combined_rainfall['DATE'], errors='coerce')
             combined_rainfall = combined_rainfall.dropna(subset=['DATE'])

        print(f"Final processed rainfall data shape: {combined_rainfall.shape}")
        print(f"Final processed tide data shape: {tide_data.shape}")

        if combined_rainfall.empty or tide_data.empty:
            print("Error: Dataframes became empty after final cleaning/date conversion.")
            return None, None

        # Report date ranges
        min_rain_date = combined_rainfall['DATE'].min()
        max_rain_date = combined_rainfall['DATE'].max()
        min_tide_date = tide_data['DATE'].min()
        max_tide_date = tide_data['DATE'].max()

        if pd.isna(min_rain_date) or pd.isna(max_rain_date): print("Warning: Could not determine rainfall date range.")
        else: print(f"Rainfall date range: {min_rain_date.date()} to {max_rain_date.date()}")

        if pd.isna(min_tide_date) or pd.isna(max_tide_date): print("Warning: Could not determine tide date range.")
        else: print(f"Tide date range: {min_tide_date.date()} to {max_tide_date.date()}")

        return tide_data, combined_rainfall

    except FileNotFoundError as e:
        print(f"Error loading data: File not found - {e}.")
        return None, None
    except Exception as e:
        print(f"Unexpected error during data loading: {e}")
        import traceback
        traceback.print_exc()
        return None, None

# Modify the process_multiple_rainfall_files function to include the export step
def process_multiple_rainfall_files(rainfall_data_paths):
    """
    Process multiple rainfall data files and combine them by priority.
    Handles different formats based on naming conventions or content sniffing.
    """
    all_rainfall_dfs = []
    print(f"\nProcessing {len(rainfall_data_paths)} rainfall data file(s)...")

    for i, file_path in enumerate(rainfall_data_paths):
        print(f"--> Processing rainfall file {i+1}/{len(rainfall_data_paths)}: '{file_path}'")
        processed_df = None
        try:
            # 1. Try to determine format based on filename convention
            filename = file_path.split('/')[-1].split('\\')[-1].lower() # Get base filename in lowercase
            if filename.startswith("daily_") and filename.endswith(".csv"):
                print("    Format detected: 'daily_*' (likely 5 columns, header)")
                processed_df = parse_daily_rainfall_file(file_path)
            elif filename.endswith("_rain.csv"):
                 print("    Format detected: '*_rain.csv' (likely 4 columns, no header)")
                 processed_df = parse_rain_suffix_file(file_path)
            else:
                # 2. If convention doesn't match, try sniffing content
                print("    Format not recognized by filename, attempting content sniffing...")
                try:
                    with open(file_path, 'r', encoding='utf-8') as f:
                        # Read a few lines to guess structure
                        lines = [f.readline().strip() for _ in range(10)]
                        # Remove empty lines
                        lines = [line for line in lines if line]
                        if not lines: raise ValueError("File seems empty.")

                        # Guess based on comma count in data lines (skip potential headers)
                        data_lines = [line for line in lines if line and not line.lower().startswith(('year', '#', '"', 'station'))] # Simple header check
                        if not data_lines: data_lines = lines # Fallback if no clear data lines found

                        comma_counts = [line.count(',') for line in data_lines]
                        if not comma_counts: raise ValueError("No commas found in sample lines.")

                        avg_commas = np.mean(comma_counts)

                        if avg_commas >= 4: # Likely 5+ columns -> daily format?
                            print("    Sniffing suggests 'daily_*' format (>=4 commas).")
                            processed_df = parse_daily_rainfall_file(file_path)
                        elif avg_commas == 3: # Likely 4 columns -> rain suffix format?
                            print("    Sniffing suggests '*_rain.csv' format (3 commas).")
                            processed_df = parse_rain_suffix_file(file_path)
                        else:
                            # print(f"    Warning: Unrecognized format based on sniffing (avg commas ~ {avg_commas:.1f}). Skipping file.")
                            # continue # Skip this file # OLD
                            pass
                except Exception as sniff_err:
                    # print(f"    Error during content sniffing for '{file_path}': {sniff_err}. Skipping file.")
                    # continue # Skip this file # OLD

                    pass
            # 3. Validate processed data
            if processed_df is not None and not processed_df.empty:
                # Basic validation: Check for DATE and Rainfall columns
                if 'DATE' in processed_df.columns and 'Rainfall' in processed_df.columns:
                    # Ensure DATE is datetime
                    if not pd.api.types.is_datetime64_any_dtype(processed_df['DATE']):
                         processed_df['DATE'] = pd.to_datetime(processed_df['DATE'], errors='coerce')
                    # Ensure Rainfall is numeric
                    if not pd.api.types.is_numeric_dtype(processed_df['Rainfall']):
                         processed_df['Rainfall'] = pd.to_numeric(processed_df['Rainfall'], errors='coerce')
                    # Drop rows where essential columns became NaT/NaN after conversion
                    processed_df.dropna(subset=['DATE', 'Rainfall'], inplace=True)

                    if not processed_df.empty:
                        processed_df['priority'] = i # Lower index = higher priority
                        all_rainfall_dfs.append(processed_df)
                        print(f"    Successfully loaded and processed. Shape: {processed_df.shape}")
                    else:
                         print("    Warning: Data became empty after essential column conversion/cleaning.")
                else:
                    print("    Warning: Processed DataFrame missing required 'DATE' or 'Rainfall' column.")
            elif processed_df is None:
                 print("    Processing function returned None.")
            else: # processed_df is empty
                print("    Warning: File processed but resulted in an empty DataFrame.")

        except FileNotFoundError:
             print(f"    Error: File not found at '{file_path}'. Skipping.")
        except Exception as e:
            print(f"    Error processing '{file_path}': {e}")
            import traceback
            traceback.print_exc() # Print detailed traceback for debugging

    # --- Combine DataFrames ---
    if not all_rainfall_dfs:
        print("\nError: No valid rainfall data could be loaded from any provided file.")
        return None

    print(f"\nCombining data from {len(all_rainfall_dfs)} successfully processed rainfall file(s)...")
    # Use concat instead of merge if they have the same columns (DATE, Rainfall, priority, etc.)
    try:
        combined_df = pd.concat(all_rainfall_dfs, ignore_index=True)
    except Exception as e:
        print(f"Error during concatenation: {e}")
        return None

    print(f"  Total rows before duplicate removal: {len(combined_df)}")
    # Sort by DATE and then priority (lower priority number means higher preference)
    combined_df.sort_values(['DATE', 'priority'], ascending=[True, True], inplace=True)

    # Keep the first occurrence of each date (which corresponds to the lowest priority number)
    rows_before_dedup = len(combined_df)
    combined_df = combined_df.drop_duplicates(subset=['DATE'], keep='first')
    rows_after_dedup = len(combined_df)
    print(f"  Removed {rows_before_dedup - rows_after_dedup} duplicate date entries based on priority.")

    # Keep track of data sources
    station_column = 'Station'
    if 'Station' not in combined_df.columns:
        combined_df[station_column] = ''
    for df in all_rainfall_dfs:
        if 'priority' in df.columns:
            priority = df['priority'].iloc[0] if not df.empty else 999
            station = extract_station_name(rainfall_data_paths[priority]) if priority < len(rainfall_data_paths) else "Unknown"
            # Mark rows in combined_df from this source
            for date in df['DATE'].unique():
                mask = (combined_df['DATE'] == date) & (combined_df['priority'] == priority)
                combined_df.loc[mask, station_column] = station

    # Drop the temporary priority column AND Completeness if it exists
    cols_to_drop = ['priority']
    if 'Completeness' in combined_df.columns:
        cols_to_drop.append('Completeness')
    combined_df = combined_df.drop(columns=cols_to_drop, errors='ignore')

    # Final check for empty dataframe
    if combined_df.empty:
        print("Error: Combined rainfall data is empty after deduplication.")
        return None

    print(f"Final combined rainfall data shape: {combined_df.shape}")
    return combined_df

def parse_rain_suffix_file(file_path):
    """
    Parse rainfall file in '*_rain.csv' format (assumed 4 cols, no header).
    Handles '******', '***', '#' indicators.
    """
    try:
        # Read with specific options for this format
        rainfall_data = pd.read_csv(
            file_path,
            encoding='utf-8',
            header=None, # Assume no header
            names=['Year', 'Month', 'Day', 'Rainfall'],
            low_memory=False
        )
        print(f"    Initial read shape: {rainfall_data.shape}")

        # Data Cleaning Steps
        # 1. Convert 'Rainfall' to string for processing markers
        rainfall_str = rainfall_data['Rainfall'].astype(str).str.strip()

        # 2. Handle missing data markers ('******', '***') -> NaN
        missing_markers = ['******', '***']
        missing_mask = rainfall_str.isin(missing_markers)
        if missing_mask.any():
            print(f"    Note: Found {missing_mask.sum()} missing data markers ({missing_markers}). Setting to NaN.")
            rainfall_data.loc[missing_mask, 'Rainfall'] = np.nan

        # 3. Handle incompleteness marker '#' - convert to NaN (MODIFIED)
        incomplete_mask = rainfall_str.str.contains('#', na=False)
        if incomplete_mask.any():
            print(f"    Note: Found {incomplete_mask.sum()} values with '#' (incomplete). Setting to NaN.")
            rainfall_data.loc[incomplete_mask, 'Rainfall'] = np.nan
            # Store completeness info
            rainfall_data['Completeness'] = 'C'  # Default complete
            rainfall_data.loc[incomplete_mask, 'Completeness'] = '#'  # Mark incomplete
            rainfall_data.loc[missing_mask, 'Completeness'] = 'M'  # Mark missing
        else:
            rainfall_data['Completeness'] = 'C'  # All seem complete if no '#' found
            rainfall_data.loc[missing_mask, 'Completeness'] = 'M'  # Mark missing

        # 4. Convert 'Rainfall' column to numeric, coercing errors
        rainfall_data['Rainfall'] = pd.to_numeric(rainfall_data['Rainfall'], errors='coerce')
        num_nan_rainfall = rainfall_data['Rainfall'].isna().sum()
        # Check if new NaNs were created beyond explicitly marked missing
        if num_nan_rainfall > missing_mask.sum() + incomplete_mask.sum():
             print(f"    Warning: {num_nan_rainfall - missing_mask.sum() - incomplete_mask.sum()} additional Rainfall values failed numeric conversion.")

        # 5. Process date columns to numeric, coercing errors
        date_cols = ['Year', 'Month', 'Day']
        for col in date_cols:
            rainfall_data[col] = pd.to_numeric(rainfall_data[col], errors='coerce')

        # 6. Drop rows with invalid dates (NaN in Year/Month/Day)
        original_rows = len(rainfall_data)
        rainfall_data.dropna(subset=date_cols, inplace=True)
        rows_dropped_date = original_rows - len(rainfall_data)
        if rows_dropped_date > 0:
            print(f"    Warning: Removed {rows_dropped_date} rows with invalid (NaN) date components.")

        # 7. Convert date columns to integer
        for col in date_cols:
            rainfall_data[col] = rainfall_data[col].astype(int)

        # 8. Create DATE column
        rainfall_data = create_date_column(rainfall_data, f"rainfall ('{file_path}')")
        if rainfall_data is None or rainfall_data.empty:
             print(f"    Error: Failed to create valid DATE column or data became empty.")
             return None

        # 9. Select and return essential columns
        essential_cols = ['DATE', 'Year', 'Month', 'Day', 'Rainfall', 'Completeness']
        # Ensure all essential columns exist before returning
        final_cols = [col for col in essential_cols if col in rainfall_data.columns]
        return rainfall_data[final_cols]

    except Exception as e:
        print(f"    Error parsing _rain.csv rainfall file '{file_path}': {e}")
        # import traceback # Uncomment for detailed debugging
        # traceback.print_exc() # Uncomment for detailed debugging
        return None
    
def parse_daily_rainfall_file(file_path):
    """
    Parse rainfall file in 'daily_*' format (assumed 5 cols, skip 3 header).
    Handles 'Trace', '******', '***', '#' indicators.
    """
    try:
        # Read with specific options for this format
        rainfall_data = pd.read_csv(
            file_path,
            encoding='utf-8',
            skiprows=3,
            names=['Year', 'Month', 'Day', 'Rainfall', 'Completeness'],
            # Use low_memory=False if dealing with mixed types or large files
            low_memory=False
        )
        print(f"    Initial read shape: {rainfall_data.shape}")

        # Data Cleaning Steps
        # 1. Convert 'Trace' to 0 in 'Rainfall' column BEFORE numeric conversion
        trace_mask = rainfall_data['Rainfall'].astype(str).str.strip() == 'Trace'
        if trace_mask.any():
            print(f"    Note: Replacing {trace_mask.sum()} 'Trace' values with 0.")
            rainfall_data.loc[trace_mask, 'Rainfall'] = 0

        # 2. Handle missing data markers ('******', '***') -> NaN
        missing_markers = ['******', '***']
        missing_mask = rainfall_data['Rainfall'].astype(str).isin(missing_markers)
        if missing_mask.any():
            print(f"    Note: Found {missing_mask.sum()} missing data markers ({missing_markers}). Setting to NaN.")
            rainfall_data.loc[missing_mask, 'Rainfall'] = np.nan
            
        # 3. Handle incompleteness marker '#' - LOGIC REMOVED / CHANGED
        #    '#' is already removed from Rainfall values. No 'Completeness' column is created.
        #    The original incomplete_mask logic is no longer needed here.

        # 4. Convert 'Rainfall' column to numeric, coercing errors
        original_nan_count_before_coerce = rainfall_data['Rainfall'].isna().sum() # NaNs from missing_markers
        rainfall_data['Rainfall'] = pd.to_numeric(rainfall_data['Rainfall'], errors='coerce')
        final_nan_count = rainfall_data['Rainfall'].isna().sum()
        
        newly_coerced_nans = final_nan_count - original_nan_count_before_coerce
        if newly_coerced_nans > 0:
             print(f"    Warning: {newly_coerced_nans} additional Rainfall values (originally not '{','.join(missing_markers)}') failed numeric conversion and were set to NaN.")

        # 5. Process date columns to numeric, coercing errors
        date_cols = ['Year', 'Month', 'Day']
        for col in date_cols:
            rainfall_data[col] = pd.to_numeric(rainfall_data[col], errors='coerce')

        # 6. Drop rows with invalid dates (NaN in Year/Month/Day)
        original_rows = len(rainfall_data)
        rainfall_data.dropna(subset=date_cols, inplace=True)
        rows_dropped_date = original_rows - len(rainfall_data)
        if rows_dropped_date > 0:
            print(f"    Warning: Removed {rows_dropped_date} rows with invalid (NaN) date components.")

        # 7. Convert date columns to integer
        for col in date_cols:
            rainfall_data[col] = rainfall_data[col].astype(int)

        # 8. Create DATE column
        rainfall_data = create_date_column(rainfall_data, f"rainfall ('{file_path}')")
        if rainfall_data is None or rainfall_data.empty:
             print(f"    Error: Failed to create valid DATE column or data became empty.")
             return None

        # 9. Select and return essential columns - 'Completeness' is NOT included
        essential_cols = ['DATE', 'Year', 'Month', 'Day', 'Rainfall']
        # Ensure all essential columns exist before returning
        final_cols = [col for col in essential_cols if col in rainfall_data.columns]
        # Check if all expected essential_cols are present
        if len(final_cols) != len(essential_cols):
            missing_essential = set(essential_cols) - set(final_cols)
            print(f"    Warning: Essential columns {missing_essential} not found in processed data for {file_path}. Returning available.")

        return rainfall_data[final_cols]

    except Exception as e:
        print(f"    Error parsing daily rainfall file '{file_path}': {e}")
        # import traceback # Uncomment for detailed debugging
        # traceback.print_exc() # Uncomment for detailed debugging
        return None

def extract_station_name(file_path):
    """
    Extract station name from file path.
    For example:
    - "QUB_tide.csv" -> "QUB"
    - "WAG_tide.csv" -> "WAG"
    - "daily_HKO_RF_ALL.csv" -> "HKO"
    """
    if not file_path:
        return "Unknown"
    
    try:
        # Get the filename without path
        filename = file_path.split('/')[-1].split('\\')[-1]  # Works for / or \ separators
        
        # Check if it's a tide file
        if "_tide" in filename.lower():
            # Extract station name from tide file (before "_tide")
            name = filename.split("_tide")[0]
            return name
            
        # Check if it's a daily rainfall file
        elif filename.startswith("daily_"):
            # Extract station name from daily file (between "daily_" and "_RF")
            if "_RF" in filename:
                name = filename.split("daily_")[1].split("_RF")[0]
                return name
            # Fallback if "_RF" not found
            else:
                name = filename.split("daily_")[1].split(".")[0]
                return name
                
        # Check if it's a rainfall file with _rain suffix
        elif "_rain" in filename.lower():
            # Extract station name from rain file (before "_rain")
            name = filename.split("_rain")[0]
            return name
            
        # Default case: just use filename without extension
        else:
            return filename.split('.')[0]
            
    except Exception as e:
        print(f"Warning: Could not extract station name from '{file_path}': {e}")
        return "Unknown"

def create_date_column(df, source_name="data"):
    """
    Create a DATE column from Year, Month, Day columns.
    Handles potential errors during conversion.
    """
    required_cols = ['Year', 'Month', 'Day']
    if not all(col in df.columns for col in required_cols):
        print(f"Error in create_date_column ({source_name}): Missing one or more required columns: {required_cols}")
        return None

    # Attempt vectorized conversion first (faster)
    try:
        # Ensure columns are suitable for datetime conversion (e.g., numeric)
        # Errors should have been coerced earlier, but double-check
        for col in required_cols:
             if not pd.api.types.is_numeric_dtype(df[col]):
                  df[col] = pd.to_numeric(df[col], errors='coerce')

        # Drop rows where date components are still NaN after coercion attempt
        original_count = len(df)
        df.dropna(subset=required_cols, inplace=True)
        if len(df) < original_count:
            print(f"Warning in create_date_column ({source_name}): Dropped {original_count - len(df)} rows with NaN in date components.")

        # Convert to integer if they are floats after coercion
        for col in required_cols:
            if pd.api.types.is_float_dtype(df[col]):
                 df[col] = df[col].astype(int)

        df['DATE'] = pd.to_datetime(df[required_cols])
        # Check if any dates failed conversion (became NaT)
        nat_count = df['DATE'].isna().sum()
        if nat_count > 0:
            print(f"Warning in create_date_column ({source_name}): {nat_count} rows resulted in NaT dates during conversion. Removing them.")
            df.dropna(subset=['DATE'], inplace=True)
        return df
    except (ValueError, TypeError, pd.errors.OutOfBoundsDatetime) as e:
        print(f"Error during vectorized date creation for {source_name}: {e}. Trying row-by-row (slower).")
        valid_indices = []
        invalid_count = 0
        for idx, row in df.iterrows():
            try:
                # Check if components are valid integers first
                year, month, day = int(row['Year']), int(row['Month']), int(row['Day'])
                # Attempt to create timestamp
                pd.Timestamp(year=year, month=month, day=day)
                valid_indices.append(idx)
            except (ValueError, TypeError, pd.errors.OutOfBoundsDatetime):
                invalid_count += 1

        if invalid_count > 0:
            print(f"Warning: Skipped {invalid_count} invalid date rows in {source_name} during row-by-row check.")

        df_filtered = df.loc[valid_indices].copy()

        if not df_filtered.empty:
            # Re-attempt conversion on the filtered data (should succeed)
            try:
                df_filtered['DATE'] = pd.to_datetime(df_filtered[['Year', 'Month', 'Day']])
                # Final check for NaT in filtered data
                if df_filtered['DATE'].isna().any():
                     print(f"Error: NaT dates found even after row-by-row validation for {source_name}. Check data.")
                     df_filtered.dropna(subset=['DATE'], inplace=True)
                return df_filtered
            except Exception as e_inner:
                 print(f"Error creating DATE column even after row-by-row filtering for {source_name}: {e_inner}")
                 return None
        else:
            print(f"Error: No valid dates found in {source_name} after filtering.")
            return None
    except Exception as e_outer:
        print(f"Unexpected error in create_date_column for {source_name}: {e_outer}")
        return None


def extract_tide_station_prefix(tide_path):
    """Extract tide station prefix from tide file path."""
    if not tide_path:
        return "Unknown"
    try:
        # Get the filename without path
        filename = tide_path.split('/')[-1].split('\\')[-1] # Works for / or \ separators
        # Extract the prefix (usually before "_tide.csv")
        parts = filename.split('_')
        if len(parts) > 1 and "tide" in parts[-1].lower():
            return parts[0]
        else:
            # Fallback: use filename without extension if no underscore found
            return filename.split('.')[0]
    except Exception:
        return "Unknown" # Fallback on any error

def export_merged_rainfall(merged_rainfall_data, tide_data_path, output_dir=""):
    """
    Export merged rainfall data to CSV file with name based on station.
    
    Parameters:
    -----------
    merged_rainfall_data : pd.DataFrame
        The merged rainfall data to export
    tide_data_path : str
        Path to the tide data file (used to extract station name)
    output_dir : str, optional
        Directory to save output file
    
    Returns:
    --------
    str
        Path to the saved file
    """
    if merged_rainfall_data is None or merged_rainfall_data.empty:
        print("Error: No rainfall data to export")
        return None
        
    # Extract station name from tide file path
    station_name = extract_station_name(tide_data_path)
    
    # Create output filename
    if output_dir and not output_dir.endswith(('/', '\\')):
        output_dir += '/'
    output_path = f"{output_dir}rainfall_all_{station_name}.csv"
    
    # Ensure DATE is in a readable format
    export_data = merged_rainfall_data.copy()
    if pd.api.types.is_datetime64_any_dtype(export_data['DATE']):
        export_data['DATE'] = export_data['DATE'].dt.strftime('%Y-%m-%d')
    
    # Sort by date before exporting
    export_data = export_data.sort_values('DATE')
    
    # Export to CSV
    export_data.to_csv(output_path, index=False, encoding='utf-8-sig')
    
    print(f"Exported merged rainfall data to '{output_path}' ({len(export_data)} records)")
    return output_path

# --- Compound Event Creation ---
def create_compound_events(tide_data, rainfall_data, rainy_day_def=0.1, rainfall_threshold_method='rainy_days'):
    """
    Merge data, calculate initial thresholds for reporting/comparison
    (but POT thresholds will be calculated later).
    Identify initial compound events based on these thresholds for clustering.
    """
    # logging omitted
    if tide_data is None or rainfall_data is None or tide_data.empty or rainfall_data.empty:
        # logging omitted
        return None, np.nan, np.nan # Return None and NaN thresholds

    # Ensure DATE columns are suitable for merging
    if not pd.api.types.is_datetime64_any_dtype(tide_data['DATE']):
        tide_data['DATE'] = pd.to_datetime(tide_data['DATE'], errors='coerce').dropna()
    if not pd.api.types.is_datetime64_any_dtype(rainfall_data['DATE']):
        rainfall_data['DATE'] = pd.to_datetime(rainfall_data['DATE'], errors='coerce').dropna()

    if tide_data.empty or rainfall_data.empty:
        # logging omitted
        return None, np.nan, np.nan

    # Merge data
    print("Merging rainfall and tide data on DATE...")
    # Select necessary columns before merge to avoid conflicts
    rain_subset = rainfall_data[['DATE', 'Rainfall', 'Year', 'Month', 'Day']].copy()
    tide_subset = tide_data[['DATE', 'Sea_Level']].copy()
    compound_data = pd.merge(rain_subset, tide_subset, on='DATE', how='inner')
    print(f"Merged data shape: {compound_data.shape}")

    if compound_data.empty:
        # logging omitted
        min_r, max_r = rainfall_data['DATE'].min(), rainfall_data['DATE'].max()
        min_t, max_t = tide_data['DATE'].min(), tide_data['DATE'].max()
        print(f"  Rainfall date range: {min_r} to {max_r}")
        print(f"  Tide date range: {min_t} to {max_t}")
        latest_start = max(min_r, min_t)
        earliest_end = min(max_r, max_t)
        print(f"  Overlap period: {latest_start} to {earliest_end}")
        if latest_start > earliest_end:
            print("  No date overlap between datasets.")
        return None, np.nan, np.nan

    # Report coverage
    date_min = compound_data['DATE'].min(); date_max = compound_data['DATE'].max()
    if pd.isna(date_min) or pd.isna(date_max):
        print("Warning: Cannot calculate data coverage due to missing dates.")
        coverage=0; total_days=0; actual_days=len(compound_data)
    else:
        total_days = (date_max - date_min).days + 1
        actual_days = len(compound_data) # Or compound_data['DATE'].nunique() for unique days
        coverage = actual_days / total_days * 100 if total_days > 0 else 0
    print(f"Merged data date range: {date_min.date()} to {date_max.date()}")
    print(f"Data coverage: {coverage:.1f}% ({actual_days} days with data within a {total_days}-day span)")

    # --- Calculate Initial Thresholds for Event Identification/Clustering ---
    # These are NOT necessarily the POT thresholds used for marginal fitting later.
    # We calculate them here to identify events for the 'cluster_compound_events' function.
    print("\nCalculating initial 95th percentile thresholds for event identification...")
    valid_sea_level = compound_data['Sea_Level'].dropna()
    valid_rainfall = compound_data['Rainfall'].dropna()

    if len(valid_sea_level) < 20 or len(valid_rainfall) < 20: # Check if enough data for percentiles
        print("Warning: Insufficient data points to reliably calculate 95th percentiles.")
        sea_level_95th = np.nan
        rainfall_95th_overall = np.nan
        rainfall_95th_rainy = np.nan
    else:
        sea_level_95th = np.percentile(valid_sea_level, 95)
        print(f"  Sea level 95th percentile (overall): {sea_level_95th:.3f} m")

        # Rainfall threshold based on method
        rainfall_95th_overall = np.percentile(valid_rainfall, 95)
        print(f"  Rainfall 95th percentile (all days): {rainfall_95th_overall:.2f} mm")

        rainy_days_data = valid_rainfall[valid_rainfall > rainy_day_def]
        if len(rainy_days_data) >= 20: # Check enough rainy days
            rainfall_95th_rainy = np.percentile(rainy_days_data, 95)
            print(f"  Rainfall 95th percentile (rainy days > {rainy_day_def} mm): {rainfall_95th_rainy:.2f} mm")
        else:
            print(f"  Warning: Insufficient rainy days ({len(rainy_days_data)}) for rainy day 95th percentile. Using overall threshold.")
            rainfall_95th_rainy = rainfall_95th_overall # Fallback

    # Choose the rainfall threshold for initial event definition based on the method
    if rainfall_threshold_method == 'rainy_days':
        initial_rainfall_threshold = rainfall_95th_rainy
        print(f"Using RAINY DAY 95th percentile for initial event definition: {initial_rainfall_threshold:.2f} mm")
    elif rainfall_threshold_method == 'all_days':
        initial_rainfall_threshold = rainfall_95th_overall
        print(f"Using ALL DAY 95th percentile for initial event definition: {initial_rainfall_threshold:.2f} mm")
    else:
        print(f"Warning: Invalid rainfall_threshold_method '{rainfall_threshold_method}'. Defaulting to 'rainy_days'.")
        initial_rainfall_threshold = rainfall_95th_rainy
        print(f"Using RAINY DAY 95th percentile for initial event definition: {initial_rainfall_threshold:.2f} mm")

    initial_sea_level_threshold = sea_level_95th
    print(f"Using Sea level 95th percentile for initial event definition: {initial_sea_level_threshold:.3f} m")

    # Handle case where thresholds couldn't be calculated
    if pd.isna(initial_rainfall_threshold) or pd.isna(initial_sea_level_threshold):
         # logging omitted
         return None, np.nan, np.nan

    # Identify initial extreme events based on these thresholds
    compound_data['Is_Extreme_Rainfall'] = compound_data['Rainfall'] >= initial_rainfall_threshold
    compound_data['Is_Extreme_Sea_Level'] = compound_data['Sea_Level'] >= initial_sea_level_threshold
    compound_data['Is_Compound_Extreme'] = (compound_data['Is_Extreme_Rainfall'] & compound_data['Is_Extreme_Sea_Level'])

    extreme_rainfall_count = compound_data['Is_Extreme_Rainfall'].sum()
    extreme_sea_level_count = compound_data['Is_Extreme_Sea_Level'].sum()
    compound_extreme_count = compound_data['Is_Compound_Extreme'].sum()
    print(f"\nInitial extreme event counts (based on 95th percentiles):")
    print(f"  Extreme rainfall events: {extreme_rainfall_count}")
    print(f"  Extreme sea level events: {extreme_sea_level_count}")
    print(f"  Initial compound extreme events (R >= {initial_rainfall_threshold:.2f} AND S >= {initial_sea_level_threshold:.3f}): {compound_extreme_count}")

    # Cluster the initial compound events to find independent representatives
    print("\nClustering initial compound extreme events to find independent representatives...")
    compound_data_clustered = cluster_compound_events(compound_data) # This adds 'Is_Rep_Compound_Extreme'

    # Return the clustered data and the initial thresholds used for clustering
    # The POT thresholds for fitting will be determined later
    return compound_data_clustered, initial_rainfall_threshold, initial_sea_level_threshold

# --- Cluster Compound Events (Keep original logic) ---
def cluster_compound_events(compound_data, min_interval=3):
    """Cluster compound floods, select most severe from each cluster."""
    data = compound_data.copy().sort_values('DATE').reset_index(drop=True)

    # Identify rows where initial compound condition is met
    compound_indices = data.index[data['Is_Compound_Extreme']].tolist()

    if not compound_indices:
        print("No initial compound events found to cluster.")
        data['Is_Rep_Compound_Extreme'] = False
        return data

    print(f"Found {len(compound_indices)} initial compound event days.")
    # Assign cluster IDs
    cluster_id = 0
    last_event_index = -min_interval # Initialize to allow the first event
    data['Cluster_ID'] = -1

    for current_idx in compound_indices:
        # Check time difference from the last assigned event in a cluster
        if current_idx - last_event_index >= min_interval:
            cluster_id += 1 # Start a new cluster
        data.loc[current_idx, 'Cluster_ID'] = cluster_id
        last_event_index = current_idx # Update the index of the last event assigned

    total_clusters = cluster_id
    print(f"Grouped into {total_clusters} initial clusters based on min_interval={min_interval} days.")

    # Select representative event for each cluster
    representative_events_indices = []
    print("Selecting representative event from each cluster based on severity score...")
    # Use overall data (excluding NaNs) for percentile calculation
    # This makes the severity score consistent across clusters
    rainfall_all_valid = data['Rainfall'].dropna()
    sealevel_all_valid = data['Sea_Level'].dropna()

    if rainfall_all_valid.empty or sealevel_all_valid.empty:
        print("Error: Cannot calculate severity scores because full dataset Rainfall or SeaLevel is empty after dropna.")
        data['Is_Rep_Compound_Extreme'] = False
        return data

    for c_id in range(1, total_clusters + 1): # Iterate through cluster IDs 1 to total_clusters
        cluster_indices = data.index[data['Cluster_ID'] == c_id].tolist()
        if not cluster_indices:
            # This shouldn't happen if indexing is correct, but check anyway
            print(f"  Warning: No events found for Cluster ID {c_id}. Skipping.")
            continue

        cluster_data = data.loc[cluster_indices].copy()

        # Calculate severity based on percentiles within the *entire* dataset
        try:
            cluster_data['Rainfall_Percentile'] = cluster_data['Rainfall'].apply(
                lambda x: stats.percentileofscore(rainfall_all_valid, x) if pd.notna(x) else 0
            )
            cluster_data['Sea_Level_Percentile'] = cluster_data['Sea_Level'].apply(
                lambda x: stats.percentileofscore(sealevel_all_valid, x) if pd.notna(x) else 0
            )
            cluster_data['Severity_Score'] = cluster_data['Rainfall_Percentile'] + cluster_data['Sea_Level_Percentile']

            # Find the index within the original dataframe 'data' corresponding to max severity
            # idxmax() called on the cluster_data slice will return an index from 'data'
            representative_index = cluster_data['Severity_Score'].idxmax()
            representative_events_indices.append(representative_index)

        except Exception as e:
            print(f"  Error calculating severity or finding max for cluster {c_id}: {e}. Skipping cluster.")
            continue

    # Mark the representative events in the original dataframe
    data['Is_Rep_Compound_Extreme'] = False
    if representative_events_indices:
        # Ensure indices are valid before using .loc
        valid_indices = [idx for idx in representative_events_indices if idx in data.index]
        if len(valid_indices) != len(representative_events_indices):
             print(f"  Warning: Some representative indices were invalid. Found {len(valid_indices)} valid.")
        if valid_indices:
            data.loc[valid_indices, 'Is_Rep_Compound_Extreme'] = True
            final_rep_count = len(valid_indices)
            print(f"Selected {final_rep_count} representative compound events after clustering.")
        else:
             print("  Warning: No valid representative events found after filtering.")
             final_rep_count = 0
    else:
        print("  No representative events selected.")
        final_rep_count = 0

    # Report counts
    original_compound_count = data['Is_Compound_Extreme'].sum()
    # final_rep_count = data['Is_Rep_Compound_Extreme'].sum() # Already calculated
    print(f"Summary: Initial compound event days = {original_compound_count}, Independent representative events = {final_rep_count}")

    return data.drop(columns=['Cluster_ID']) # Remove temporary cluster ID column


# Calculate probability based on distribution
def calculate_probability(value, dist_info):
    """Internal helper."""
    dist_name = dist_info['distribution']
    params = dist_info['params']

    try:
        if dist_name == 'genpareto':
            threshold = dist_info.get('threshold', 0)
            if value <= threshold:
                return 0.0
            exceedance = value - threshold
            p_excess = stats.genpareto.cdf(exceedance, *params)
            p_exceed_threshold = len(dist_info.get('exceedance_values_for_calc_prob', [])) / dist_info.get('num_samples_for_pot_calc_prob', 1000)
            return 1.0 - p_exceed_threshold * (1.0 - p_excess)

        return getattr(stats, dist_name).cdf(value, *params)
    except Exception:
        return 0.5

# Probability integral transform
def probability_transform(data, fits, compound_data):
    """
    Transform data to uniform [0,1] interval using rank method
    """
    u_data = pd.DataFrame(index=data.index)
    
    for column in fits.keys():
        # Use rank method to calculate cumulative probability
        values = data[column].values
        ranks = stats.rankdata(values)
        # Use Weibull estimator: (rank - 0.5) / n
        # This ensures values are uniformly distributed in (0,1)
        n = len(values)
        u_values = (ranks - 0.5) / n
        
        # Ensure strictly in (0,1) interval, avoid 0 and 1
        u_values = np.clip(u_values, 0.001, 0.999)
        
        u_data[column] = u_values
    
    # Diagnostic output
    for column in fits.keys():
        min_val = u_data[column].min()
        max_val = u_data[column].max()
        print(f"Transformed {column} range: [{min_val:.6f}, {max_val:.6f}]")
    
    return u_data

# Perform Kolmogorov-Smirnov test
def perform_kolmogorov_smirnov_test(data, copula_model):
    """
    Perform Kolmogorov-Smirnov test to assess goodness of fit
    """
    n = len(data)
    u = data[:, 0]
    v = data[:, 1]
    
    # Calculate empirical Copula CDF
    emp_cop = np.zeros(n)
    for i in range(n):
        emp_cop[i] = np.mean((u <= u[i]) & (v <= v[i]))
    
    # Calculate theoretical Copula CDF
    if copula_model is None:
        # Independent Copula
        theo_cop = u * v
    else:
        # Use copulae library's CDF method
        theo_cop = np.zeros(n)
        for i in range(n):
            try:
                theo_cop[i] = copula_model.cdf(data[i])
            except Exception as e:
                # Fall back to independent Copula if CDF calculation fails
                theo_cop[i] = u[i] * v[i]
    
    # Calculate Kolmogorov-Smirnov statistic (maximum absolute difference)
    ks_stat = np.max(np.abs(emp_cop - theo_cop))
    
    # Calculate p-value based on empirical formula
    p_value = np.exp(-2 * n * ks_stat ** 2)
    # p_value = min(max(p_value, 0.001), 0.999)  # Limit to reasonable range # REMOVED
    
    return p_value, ks_stat

# Perform Cramer-von Mises test
def perform_cramer_von_mises_test(data, copula_model):
    """
    Perform Cramer-von Mises test to assess goodness of fit
    """
    n = len(data)
    u = data[:, 0]
    v = data[:, 1]
    
    # Calculate empirical Copula
    emp_cop = np.zeros(n)
    for i in range(n):
        emp_cop[i] = np.mean((u <= u[i]) & (v <= v[i]))
    
    # Calculate theoretical Copula
    if copula_model is None:
        # Independent Copula
        theo_cop = u * v
    else:
        # Use copulae library's CDF method
        theo_cop = np.zeros(n)
        for i in range(n):
            try:
                theo_cop[i] = copula_model.cdf(data[i])
            except Exception as e:
                # Fall back to independent Copula if CDF calculation fails
                theo_cop[i] = u[i] * v[i]
    
    # Calculate Cramer-von Mises statistic
    cm_stat = np.sum((emp_cop - theo_cop) ** 2)
    
    # Removed p-value calculation for CM test as it was simplified and potentially problematic.
    # The statistic itself can be used for comparison or more complex p-value lookup if needed.
    # If a p-value is strictly required here, a more robust statistical method or library 
    # for its calculation based on the CM statistic and sample size would be needed.
    # For now, we return NaN for p-value to indicate it's not calculated by this simplified version.
    p_value = np.nan # Placeholder, as robust p-value calculation is non-trivial
    
    return p_value, cm_stat


# Calculate return period
def calculate_return_period(compound_data, rainfall_col, sea_level_col, copula_results, best_copula, fits):
    """
    Calculate compound extreme event return periods using theoretical method
    T(R > r, S > s) = 1 / P(R > r, S > s)
    where P(R > r, S > s) = 1 - F_R(r) - F_S(s) + F(r,s)
    """
    # Select compound extreme events
    compound_events = compound_data[compound_data['Is_Rep_Compound_Extreme'] == True].copy()
    print(f"Analyzing {len(compound_events)} compound extreme events for return periods")
    
    copula_model = copula_results[best_copula]['model']
    # Calculate data coverage in years (for information only)
    n_years = (compound_data['DATE'].max() - compound_data['DATE'].min()).days / 365.25
    print(f"Data covers {n_years:.1f} years")
    
    results = []
    
    for idx, row in compound_events.iterrows():
        r = row[rainfall_col]
        s = row[sea_level_col]
        month = row['Month']
        year = row['Year']
        
        # Calculate rainfall marginal CDF
        rainfall_fit = fits[rainfall_col]
        rainfall_cdf = calculate_probability(r, rainfall_fit)
        
        # Calculate sea level marginal CDF
        sea_level_fit = fits[sea_level_col]
        sea_level_cdf = calculate_probability(s, sea_level_fit)
        
        # Calculate joint CDF
        # Ensure marginal CDFs are strictly in (0,1) for copula calculations
        u_rain_cdf = np.clip(rainfall_cdf, 0.001, 0.999)
        v_sl_cdf = np.clip(sea_level_cdf, 0.001, 0.999)

        if copula_model is None: # Independence case
            joint_cdf_val = u_rain_cdf * v_sl_cdf
        else:
            # Check if the best copula is Survival Clayton fitted via manual transform
            is_survival_clayton_manual = (best_copula == 'survival_clayton' and 
                                          copula_results[best_copula].get('method') == 'manual_transform')

            if is_survival_clayton_manual:
                # For C_sc(u,v) = u + v - 1 + C_std(1-u, 1-v)
                u_transformed = 1 - u_rain_cdf
                v_transformed = 1 - v_sl_cdf
                # copula_model is C_std in this case
                c_std_val = copula_model.cdf([u_transformed, v_transformed]) 
                joint_cdf_val = u_rain_cdf + v_sl_cdf - 1 + c_std_val
            else:
                # For all other standard copulas (or directly supported survival copulas if any)
                try:
                    joint_cdf_val = copula_model.cdf([u_rain_cdf, v_sl_cdf])
                except Exception as e:
                    print(f"Error in Copula CDF calculation for event {row['DATE']}: {e}. Defaulting to independence.")
                    joint_cdf_val = u_rain_cdf * v_sl_cdf # Fallback to independence
        
        # Calculate exceedance probability: P(R > r, S > s) = 1 - F_R(r) - F_S(s) + F(r,s)
        # Note: F_R(r) is rainfall_cdf (u_rain_cdf), F_S(s) is sea_level_cdf (v_sl_cdf), F(r,s) is joint_cdf_val
        exceedance_prob = 1 - u_rain_cdf - v_sl_cdf + joint_cdf_val
        
        # Handle potential numerical errors, ensure probability in reasonable range
        exceedance_prob = max(exceedance_prob, 1e-6)
        exceedance_prob = min(exceedance_prob, 0.999)
        
        # Calculate theoretical return period (years) = 1 / exceedance probability
        return_period = 1 / exceedance_prob
        
        results.append({
            'Date': row['DATE'],
            'Year': year,
            'Month': month,
            'Rainfall': r,
            'Sea_Level': s,
            'Rainfall_CDF': rainfall_cdf,
            'Sea_Level_CDF': sea_level_cdf,
            'Joint_CDF': joint_cdf_val,
            'Exceedance_Prob': exceedance_prob,
            'Return_Period': return_period,
            'Is_Compound_Extreme': True
        })
    
    return pd.DataFrame(results)

def calculate_contours(compound_data, rainfall_col, sea_level_col, copula_results, best_copula, fits, 
                      rainfall_threshold, sea_level_threshold, 
                      x_min, x_max, y_min, y_max, npoints=200):
    """Internal helper."""
    # logging omitted
    
    # Create grid based on the specified ranges
    s = np.linspace(x_min, x_max, npoints)
    r = np.linspace(y_min, y_max, npoints)
    S, R = np.meshgrid(s, r)
    
    # Matrix to store joint return periods
    Z = np.zeros_like(R)
    
    # Get best Copula model
    if best_copula in copula_results:
        copula_model = copula_results[best_copula]['model']
    else:
        print(f"Warning: Copula model '{best_copula}' not found, using independent Copula")
        copula_model = None
    
    # Calculate theoretical return period for each grid point
    for i in range(len(r)):
        for j in range(len(s)):
            r_val = R[i, j]
            s_val = S[i, j]
            
            # Calculate rainfall marginal CDF
            rainfall_fit = fits[rainfall_col]
            rainfall_cdf = calculate_probability(r_val, rainfall_fit)
            
            # Calculate sea level marginal CDF
            sea_level_fit = fits[sea_level_col]
            sea_level_cdf = calculate_probability(s_val, sea_level_fit)
            
            # Calculate joint CDF
            # Ensure marginal CDFs are strictly in (0,1) for copula calculations
            u_rain_cdf_pt = np.clip(rainfall_cdf, 0.001, 0.999)
            v_sl_cdf_pt = np.clip(sea_level_cdf, 0.001, 0.999)

            if copula_model is None: # Independence case
                joint_cdf_val_pt = u_rain_cdf_pt * v_sl_cdf_pt
            else:
                is_survival_clayton_manual_pt = (best_copula == 'survival_clayton' and 
                                               copula_results[best_copula].get('method') == 'manual_transform')
                
                if is_survival_clayton_manual_pt:
                    u_transformed_pt = 1 - u_rain_cdf_pt
                    v_transformed_pt = 1 - v_sl_cdf_pt
                    c_std_val_pt = copula_model.cdf([u_transformed_pt, v_transformed_pt])
                    joint_cdf_val_pt = u_rain_cdf_pt + v_sl_cdf_pt - 1 + c_std_val_pt
                else:
                    try:
                        joint_cdf_val_pt = copula_model.cdf([u_rain_cdf_pt, v_sl_cdf_pt])
                    except:
                        joint_cdf_val_pt = u_rain_cdf_pt * v_sl_cdf_pt # Fallback
            
            # Calculate exceedance probability: P(R > r, S > s) = 1 - F_R(r) - F_S(s) + F(r,s)
            exceedance_prob = 1 - u_rain_cdf_pt - v_sl_cdf_pt + joint_cdf_val_pt
            
            # Handle potential numerical errors, ensure probability in reasonable range
            exceedance_prob = max(exceedance_prob, 1e-10)
            
            # Calculate theoretical return period (years) = 1 / exceedance probability
            return_period = 1 / exceedance_prob
            
            # Store return period in Z matrix
            Z[i, j] = return_period
    
    # Smooth Z to reduce boundary effects
    from scipy.ndimage import gaussian_filter
    Z_smooth = gaussian_filter(Z, sigma=0.5)  # Light smoothing
    
    print(f"Z matrix range: {Z_smooth.min():.2f} to {Z_smooth.max():.2f}")
    
    return R, S, Z_smooth

def find_maximum_density_points(R, S, Z, copula_results, best_copula, fits, rainfall_col, sea_level_col, 
                                return_periods=[5, 10, 20, 50, 100, 200, 500, 1000]):
    """Internal helper."""
    import numpy as np
    import pandas as pd
    from scipy import interpolate
    from scipy import stats
    
    copula_model = copula_results[best_copula]['model']
    
    contour_func = interpolate.LinearNDInterpolator(
        points=(S.flatten(), R.flatten()),
        values=Z.flatten()
    )
    
    max_density_points = []
    
    def calculate_marginal_pdf(value, dist_info):
        """Internal helper."""
        dist_name = dist_info['distribution']
        params = dist_info['params']
        try:
            if dist_name == 'genpareto':
                threshold = dist_info.get('threshold', 0)
                if value <= threshold:
                    return 0.0
                return stats.genpareto.pdf(value - threshold, *params)
            return getattr(stats, dist_name).pdf(value, *params)
        except Exception:
            return 1.0
    
    for T in return_periods:
        # logging omitted
        
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(10, 8))
        
        Z_numeric = np.array(Z, dtype=float)
        
        cs = ax.contour(S, R, Z_numeric, levels=[T])
        plt.close(fig)
        
        contour_paths = cs.collections[0].get_paths()
        if not contour_paths:
            # logging omitted
            continue
            
        all_vertices = []
        for path in contour_paths:
            vertices = path.vertices
            sea_levels = vertices[:, 0]
            rainfalls = vertices[:, 1]
            points = np.column_stack((sea_levels, rainfalls))
            all_vertices.append(points)
        
        if not all_vertices:
            # logging omitted
            continue
            
        all_points = np.vstack(all_vertices)
        
        pdf_values = np.zeros(len(all_points))
        
        for i, (s_val, r_val) in enumerate(all_points):
            rainfall_cdf = calculate_probability(r_val, fits[rainfall_col])
            sea_level_cdf = calculate_probability(s_val, fits[sea_level_col])
            
            u = np.clip(rainfall_cdf, 0.001, 0.999)
            v = np.clip(sea_level_cdf, 0.001, 0.999)
            
            if copula_model is None:
                rainfall_pdf = calculate_marginal_pdf(r_val, fits[rainfall_col])
                sea_level_pdf = calculate_marginal_pdf(s_val, fits[sea_level_col])
                pdf_values[i] = rainfall_pdf * sea_level_pdf
            else:
                try:
                    is_survival_clayton = (best_copula == 'survival_clayton' and 
                                           copula_results[best_copula].get('method') == 'manual_transform')
                    
                    if is_survival_clayton:
                        u_transformed = 1 - u
                        v_transformed = 1 - v
                        copula_pdf = copula_model.pdf([u_transformed, v_transformed])
                    else:
                        copula_pdf = copula_model.pdf([u, v])
                        
                    rainfall_pdf = calculate_marginal_pdf(r_val, fits[rainfall_col])
                    sea_level_pdf = calculate_marginal_pdf(s_val, fits[sea_level_col])
                    
                    pdf_values[i] = copula_pdf * rainfall_pdf * sea_level_pdf
                except Exception as e:
                    # logging omitted
                    pdf_values[i] = 0
        
        if np.all(pdf_values == 0):
            # logging omitted
            continue
            
        max_idx = np.argmax(pdf_values)
        max_point = all_points[max_idx]
        max_density_sea_level, max_density_rainfall = max_point
        
        point_return_period = contour_func(max_density_sea_level, max_density_rainfall)
        
        max_density_points.append({
            'ReturnPeriod': T,
            'SeaLevel_m': max_density_sea_level,
            'Rainfall_mm': max_density_rainfall,
            'PDF_Value': pdf_values[max_idx],
            'Verified_ReturnPeriod': float(point_return_period)
        })
        
        # logging omitted
    
    if max_density_points:
        results_df = pd.DataFrame(max_density_points)
        results_df.to_csv('maximum_density_points.csv', index=False, encoding='utf-8-sig')
        # logging omitted
        return results_df
    else:
        # logging omitted
        return pd.DataFrame()

def get_storm_names():
    """Return mapping of dates to storm names"""
    return {
        '1962-09-01': 'Wanda', 
        '2018-09-16': 'Mangkhut',
        '1989-05-20': 'Brenda',
        '2017-08-23': 'Hato',
        '1976-08-24': 'Ellen',
        '2001-07-06': 'Utor',
        '1960-06-09': 'Mary',
        '2008-09-24': 'Hagupit',
        '1995-08-12': 'Helen',
        '2009-09-15': 'Koppu',
        '1979-08-02': 'Hope'
    }

def plot_return_periods(ax, R, S, Z, compound_data, return_periods, 
                        rainfall_threshold, sea_level_threshold, 
                        rainfall_col='Rainfall', sea_level_col='Sea_Level'):
    """
    Modified return periods contour visualization function with a simple legend in the top left corner
    """
    # Define uniform color for points and background
    uniform_color = (80/255, 29/255, 138/255)  # Deep purple
    contour_color = '#46717C'  # Blue-green color for contours
    
    # Add background color matching the points but with high transparency
    ax.add_patch(plt.Rectangle((sea_level_threshold, rainfall_threshold), 
                              5.0-sea_level_threshold, 600-rainfall_threshold, 
                              facecolor=uniform_color, alpha=0.15, zorder=0))
    
    # Set contour levels
    levels = [5, 10, 20, 50, 100, 200, 500, 1000]
    
    # Draw contours - with swapped axes and modified colors/thickness
    cs = ax.contour(S, R, Z, levels=levels, colors=contour_color, linewidths=1.8, 
                   corner_mask=True, antialiased=True)
    
    # Add contour labels with adjusted positions
    labels = ax.clabel(cs, inline=False, fontsize=16, fmt='%d', colors=contour_color)
    
    # Adjust the position of all labels, offset upward
    for label in labels:
        pos = label.get_position()
        label_text = label.get_text()
        
        # Special handling for labels with specific values
        if label_text == '5':
            label.set_position((pos[0] + 0.01, pos[1] + 10))
        elif label_text == '10':
            label.set_position((pos[0] + 0.01, pos[1] + 12))
        elif label_text == '20':
            label.set_position((pos[0] - 0.01, pos[1] + 15))
        elif label_text == '50':
            label.set_position((pos[0] + 0.02, pos[1] + 8))
        elif label_text == '100':
            label.set_position((pos[0] + 0.02, pos[1] + 10))
        elif label_text == '200':
            label.set_position((pos[0] + 0.01, pos[1] + 12))
        elif label_text == '500':
            label.set_position((pos[0] + 0.01, pos[1] + 8))
        elif label_text == '1000':
            label.set_position((pos[0] + 0.02, pos[1] + 10))
        else:
            # Normal offset for other labels
            label.set_position((pos[0], pos[1] + 10))
    
    # Draw background points and extreme event points
    background_points = compound_data[~compound_data['Is_Compound_Extreme']]
    ax.scatter(
        background_points[sea_level_col],
        background_points[rainfall_col],
        c='lightgray',
        s=20,
        alpha=0.6,              # Reduced opacity
        edgecolors='gray',      # Add border
        linewidths=0.5,         # Border width
    )
    
    # Select compound extreme events
    compound_extreme = compound_data[compound_data['Is_Rep_Compound_Extreme'] == True]
    
    # Plot compound extreme events with uniform color and border
    scatter = ax.scatter(
        compound_extreme[sea_level_col],
        compound_extreme[rainfall_col],
        c=[uniform_color],
        s=50,
        alpha=0.6,              # Slight transparency
        edgecolors=None,        # No border
        linewidths=1.0,         # Border width
        zorder=10,
    )
    
    # Create a simple legend for return period contours only
    import matplotlib.lines as mlines
    
    # Single line for the legend
    contour_line = mlines.Line2D([], [], color=contour_color, linewidth=1.8, 
                              label='Return Period')
    
    # Add the simple legend to the top left corner
    #ax.legend(handles=[contour_line], loc='upper right', 
    #          fontsize=16, framealpha=0.8)
    
    # Get storm names
    storm_names = get_storm_names()
    
    # Label specific storms of interest with improved visibility
    custom_offsets = {
        'Wanda': (-40, 20),
        'Mangkhut': (-45, 18),
        'Brenda': (25, -25),
        'Hato': (10, 10),
        'Ellen': (15, 5),
        'Utor': (15, 10),
        'Mary': (15, 0),
        'Hagupit': (10, 5),
        'Helen': (10, 40),
        'Koppu': (29, -25),
        'Hope': (25,-10)
    }
    
    # Light pink color for storm labels
    storm_label_color = (238/255, 140/255, 125/255)
    
    for date_str, name in storm_names.items():
        date = pd.to_datetime(date_str)
        matching_data = compound_data[compound_data['DATE'] == date]
        
        if len(matching_data) > 0:
            row = matching_data.iloc[0]
            if row['Is_Compound_Extreme']:  # Check if it's a compound extreme event
                # Get custom offset for this storm or use default
                x_offset, y_offset = custom_offsets.get(name, (15, 5))
                
                ax.annotate(
                    name,
                    (row[sea_level_col], row[rainfall_col]),
                    xytext=(x_offset, y_offset),
                    textcoords='offset points',
                    fontsize=13,              # Increase font size
                    color=storm_label_color,  # Light pink color
                    fontweight='bold',        # Bold font
                    bbox=dict(boxstyle="round,pad=0.3", fc="white", ec=storm_label_color, alpha=0.7),
                    arrowprops=dict(
                        arrowstyle='-',       # Simple line, no arrow
                        color=storm_label_color,  # Use same color as label
                        lw=1.5,               # Line width
                        connectionstyle='arc3,rad=0'  # Straight line connection
                    ),
                    zorder=12  # Ensure labels are drawn on top
                )
    
    # Set axis ranges
    ax.set_xlim(sea_level_threshold, 5.0)
    ax.set_ylim(rainfall_threshold, 500)
       
    # Set x-axis ticks to one decimal place
    import matplotlib.ticker as ticker
    ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%.1f'))
    
    # Set axis labels and title with bold font
    ax.set_xlabel('Maximum Water Level $($m$)$', fontsize=20, family='Arial')
    ax.set_ylabel('24 Hours of Cumulative Rainfall $($mm$)$', fontsize=20, family='Arial')
    
    # Add grid
    ax.grid(True, linestyle='--', alpha=0.2, color='gray')
    
    # Add subplot label (b) to top right corner
    #ax.text(0.98, 0.98, '(b)', transform=ax.transAxes, 
    #        fontsize=24, fontweight='bold', ha='right', va='top')
    
    # Make the borders (spines) thicker
    for spine in ax.spines.values():
        spine.set_linewidth(1.8)
    
    # Adjust tick width, length, and label size
    ax.tick_params(axis='both', which='major', width=1.8, length=6, labelsize=20)

def plot_maximum_density_points(ax, max_density_df, color='red', marker='*', size=150):
    """Internal helper."""
    if max_density_df.empty:
        # logging omitted
        return
        
    scatter = ax.scatter(
        max_density_df['SeaLevel_m'],
        max_density_df['Rainfall_mm'],
        c=color,
        marker=marker,
        s=size,
        alpha=0.8,
        edgecolors='black',
        linewidths=1.0,
        zorder=20
    )
    
    ax.legend(loc='upper right', fontsize=12)

def export_specific_contours(R, S, Z, return_periods=[5, 10, 20, 50, 100, 200, 500, 1000], filename="specific_contours.csv"):
    """
    Extract and export contour lines for specific return periods.
    
    Parameters:
    R: 2D array - Rainfall values from meshgrid
    S: 2D array - Sea level values from meshgrid
    Z: 2D array - Return period values
    return_periods: list - List of specific return period values to extract
    filename: str - Name of the output file
    
    Returns:
    str - Path to the saved file
    """
    import pandas as pd
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.contour import QuadContourSet
    from matplotlib import cm
    
    # Create a figure (but don't display it)
    fig, ax = plt.subplots(figsize=(10, 8))
    
    # Generate contours for specified levels
    contour_set = ax.contour(S, R, Z.T, levels=return_periods)
    
    # Close the figure since we only need the contour data
    plt.close(fig)
    
    # Extract contour data
    all_contour_data = []
    
    for i, return_period in enumerate(return_periods):
        # Get contour paths for this level
        paths = contour_set.collections[i].get_paths()
        
        for path in paths:
            vertices = path.vertices
            # Extract x (sea level) and y (rainfall) coordinates
            sea_levels = vertices[:, 0]
            rainfalls = vertices[:, 1]
            
            # Create data for this contour
            contour_points = pd.DataFrame({
                'Rainfall_mm': rainfalls,
                'SeaLevel_m': sea_levels,
                'ReturnPeriod_years': return_period
            })
            
            all_contour_data.append(contour_points)
    
    # Combine all contour data
    if all_contour_data:
        contour_df = pd.concat(all_contour_data, ignore_index=True)
        
        # Round values for better readability
        contour_df = contour_df.round({
            'Rainfall_mm': 2,
            'SeaLevel_m': 3,
            'ReturnPeriod_years': 0
        })
        
        # Export to CSV
        contour_df.to_csv(filename, index=False, encoding='utf-8-sig')
        print(f"Specific contour lines exported to {filename} ({len(contour_df)} points across {len(return_periods)} contour lines)")
        return filename
    else:
        print("No contour data found for the specified return periods.")
        return None

def calculate_contour_axis_intersections(R, S, Z, return_periods=[5, 10, 20, 50, 100, 200, 500, 1000]):
    """Internal helper."""
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    from scipy import interpolate
    
    intersections_data = []
    
    # logging omitted
    
    for T in return_periods:
        # logging omitted
        
        fig, ax = plt.subplots(figsize=(8, 6))
        
        try:
            cs = ax.contour(S, R, Z, levels=[T])
            
            x_axis_intersections = []
            y_axis_intersections = []
            
            all_vertices = []
            
            try:
                if hasattr(cs, 'collections') and len(cs.collections) > 0:
                    contour_paths = cs.collections[0].get_paths()
                    for path in contour_paths:
                        vertices = path.vertices
                        if len(vertices) > 0:
                            all_vertices.append(vertices)
                else:
                    raise AttributeError("collections not available")
            except (AttributeError, IndexError):
                try:
                    if hasattr(cs, 'allsegs') and len(cs.allsegs) > 0:
                        level_segments = cs.allsegs[0]
                        for segment in level_segments:
                            if len(segment) > 0:
                                vertices = np.array(segment)
                                all_vertices.append(vertices)
                    else:
                        raise AttributeError("allsegs not available")
                except (AttributeError, IndexError):
                    try:
                        for collection in cs.collections:
                            for path in collection.get_paths():
                                vertices = path.vertices
                                if len(vertices) > 0:
                                    all_vertices.append(vertices)
                    except Exception as e3:
                        # logging omitted
                        plt.close(fig)
                        continue
            
            plt.close(fig)
            
            if not all_vertices:
                # logging omitted
                continue
            
            for vertices in all_vertices:
                sea_levels = vertices[:, 0]
                rainfalls = vertices[:, 1]
                
                min_rainfall_idx = np.argmin(rainfalls)
                min_rainfall = rainfalls[min_rainfall_idx]
                
                if len(rainfalls) > 1:
                    rainfall_range = np.max(rainfalls) - np.min(rainfalls)
                    if rainfall_range > 5:
                        if min_rainfall < (np.min(rainfalls) + rainfall_range * 0.25):
                            x_intersection = sea_levels[min_rainfall_idx]
                            if x_intersection > 0:
                                x_axis_intersections.append(x_intersection)
                
                min_sealevel_idx = np.argmin(sea_levels)
                min_sealevel = sea_levels[min_sealevel_idx]
                
                if len(sea_levels) > 1:
                    sealevel_range = np.max(sea_levels) - np.min(sea_levels)
                    if sealevel_range > 0.05:
                        if min_sealevel < (np.min(sea_levels) + sealevel_range * 0.30):
                            y_intersection = rainfalls[min_sealevel_idx]
                            if y_intersection > 0:
                                y_axis_intersections.append(y_intersection)
            
            x_intersection_final = min(x_axis_intersections) if x_axis_intersections else None
            y_intersection_final = min(y_axis_intersections) if y_axis_intersections else None
            
            if x_intersection_final is None or y_intersection_final is None:
                for vertices in all_vertices:
                    sea_levels = vertices[:, 0]
                    rainfalls = vertices[:, 1]
                    
                    if x_intersection_final is None and len(rainfalls) > 0:
                        min_rainfall_idx = np.argmin(rainfalls)
                        x_backup = sea_levels[min_rainfall_idx]
                        if x_backup > 0:
                            x_intersection_final = x_backup
                    
                    if y_intersection_final is None and len(sea_levels) > 0:
                        min_sealevel_idx = np.argmin(sea_levels)
                        y_backup = rainfalls[min_sealevel_idx]
                        if y_backup > 0:
                            y_intersection_final = y_backup
            
            notes = []
            if x_intersection_final is None:
                notes.append("No X-axis intersection found")
            elif len(x_axis_intersections) == 0:
                notes.append("X-axis intersection found using backup method")
            
            if y_intersection_final is None:
                notes.append("No Y-axis intersection found")
            elif len(y_axis_intersections) == 0:
                notes.append("Y-axis intersection found using backup method")
            
            intersections_data.append({
                'ReturnPeriod_years': T,
                'SeaLevel_at_zero_rainfall_m': x_intersection_final,
                'Rainfall_at_zero_sealevel_mm': y_intersection_final,
                'Notes': '; '.join(notes) if notes else 'Both intersections found'
            })
            
            # logging omitted
            
        except Exception as e:
            # logging omitted
            plt.close(fig)
            continue
    
    if intersections_data:
        results_df = pd.DataFrame(intersections_data)
        # logging omitted
        return results_df
    else:
        # logging omitted
        return pd.DataFrame()

# Export event results
def export_events_table(return_periods):
    """
    Export all compound event calculation results to CSV file
    """
    # Sort results by return period in descending order
    sorted_results = return_periods.sort_values('Return_Period', ascending=False)
    
    # Select relevant columns for output
    output_columns = ['Date', 'Year', 'Month', 'Rainfall', 'Sea_Level', 'Exceedance_Prob', 'Return_Period']
    output_df = sorted_results[output_columns]
    
    # Format date column for better readability
    output_df['Date'] = output_df['Date'].dt.strftime('%Y-%m-%d')
    
    # Round numeric columns for better readability
    output_df['Exceedance_Prob'] = output_df['Exceedance_Prob'].round(6)
    output_df['Return_Period'] = output_df['Return_Period'].round(2)
    
    # Export to CSV
    output_df.to_csv('compound_events_results.csv', index=False)
    print(f"Exported {len(output_df)} compound events to 'compound_events_results.csv'")
    
    # Return the formatted dataframe for display
    return output_df

def create_multi_station_figure(station_configs, figsize=(18, 12)):
    """Internal helper."""
    # Set the default font to Arial
    plt.rcParams['font.family'] = 'Arial'
    
    fig, axes = plt.subplots(2, 3, figsize=figsize, dpi=600)
    axes = axes.flatten()
    
    for idx, config in enumerate(station_configs):
        if idx >= 6:
            break
            
        ax = axes[idx]
        station_name = config['name']
        tide_file = config['tide_file']
        rainfall_files = config['rainfall_files']
        title = config.get('title', f'Station {station_name}')
        
        x_max = config.get('x_max', 5.0)
        y_max = config.get('y_max', 600)
        
        # logging omitted
        # logging omitted
        
        try:
            offset_file = config.get('offset_file', None)
            if offset_file:
                # logging omitted
                pass
            else:
                # logging omitted
            
                pass
            tide_data, rainfall_data = load_data(tide_file, rainfall_files, export_rainfall=False, 
                                               offset_file_path=offset_file)
            
            if tide_data is not None and offset_file and 'Applied_Offset' in tide_data.columns:
                applied_corrections = (tide_data['Applied_Offset'] != 0).sum()
                if applied_corrections > 0:
                    export_corrected_tide_data(tide_data, station_name)
                    # logging omitted
            
            if tide_data is None or rainfall_data is None:
                # logging omitted
                ax.text(0.5, 0.5, f'Station {station_name}\nData Loading Failed', 
                       ha='center', va='center', transform=ax.transAxes, fontsize=16)
                ax.set_title(title, fontsize=18, fontweight='bold')
                continue
            
            compound_data, rainfall_threshold, sea_level_threshold = create_compound_events(
                tide_data, rainfall_data, 
                rainy_day_def=0.1,
                rainfall_threshold_method='all_days'
            )
            
            if compound_data is None:
                # logging omitted
                ax.text(0.5, 0.5, f'Station {station_name}\nCompound Events Creation Failed', 
                       ha='center', va='center', transform=ax.transAxes, fontsize=16)
                ax.set_title(title, fontsize=18, fontweight='bold')
                continue
            
            copula_results, best_copula, fits = fit_copula_models_selective(
                compound_data, 'Rainfall', 'Sea_Level', 
                rainfall_threshold=rainfall_threshold, 
                sea_level_threshold=sea_level_threshold,
                excluded_distributions=config.get('excluded_distributions'),
            )
            
            R, S, Z = calculate_contours(
                compound_data, 'Rainfall', 'Sea_Level',
                copula_results, best_copula, fits,
                rainfall_threshold=rainfall_threshold,
                sea_level_threshold=sea_level_threshold,
                x_min=sea_level_threshold, x_max=x_max,
                y_min=rainfall_threshold, y_max=y_max,
                npoints=200
            )
            
            max_density_df = find_maximum_density_points(
                R, S, Z, copula_results, best_copula, fits,
                'Rainfall', 'Sea_Level'
            )
            
            if max_density_df is not None and not max_density_df.empty:
                csv_filename = f'maximum_density_points_{station_name}.csv'
                max_density_df.to_csv(csv_filename, index=False, encoding='utf-8-sig')
                # logging omitted
            else:
                # logging omitted
            
                pass
            try:
                intersections_df = calculate_contour_axis_intersections(
                    R, S, Z, 
                    return_periods=[5, 10, 20, 50, 100, 200, 500, 1000]
                )
                
                if intersections_df is not None and not intersections_df.empty:
                    intersections_df['Station'] = station_name
                    intersections_df = intersections_df[['Station', 'ReturnPeriod_years', 'SeaLevel_at_zero_rainfall_m', 'Rainfall_at_zero_sealevel_mm', 'Notes']]
                    
                    intersections_filename = f'contour_axis_intersections_{station_name}.csv'
                    intersections_df.to_csv(intersections_filename, index=False, encoding='utf-8-sig')
                    # logging omitted
                else:
                    # logging omitted
                    
                    pass
            except Exception as intersection_error:
                # logging omitted
            
                pass
            plot_return_periods_single_station(
                ax, R, S, Z, compound_data, max_density_df,
                rainfall_threshold, sea_level_threshold,
                station_name, title,
                x_max=x_max, y_max=y_max
            )
            
        except Exception as e:
            # logging omitted
            ax.text(0.5, 0.5, f'Station {station_name}\nProcessing Error', 
                   ha='center', va='center', transform=ax.transAxes, fontsize=16)
            ax.set_title(title, fontsize=18, fontweight='bold')
    
    for idx in range(len(station_configs), 6):
        axes[idx].set_visible(False)
    
    plt.subplots_adjust(
        left=0.08,
        right=0.95,
        top=0.92,
        bottom=0.08,
        wspace=0.35,
        hspace=0.45
    )
    
    fig.savefig('multi_station_return_periods.tiff', format='tiff', dpi=600, bbox_inches='tight')
    fig.savefig('multi_station_return_periods.png', dpi=600, bbox_inches='tight')
    # logging omitted
    
    return fig

def find_maximum_density_points(R, S, Z, copula_results, best_copula, fits, rainfall_col, sea_level_col, 
                                return_periods=[5, 10, 20, 50, 100, 200, 500, 1000]):
    """Internal helper."""
    import numpy as np
    import pandas as pd
    from scipy import interpolate
    from scipy import stats
    import matplotlib.pyplot as plt
    
    copula_model = copula_results[best_copula]['model']
    
    contour_func = interpolate.LinearNDInterpolator(
        points=(S.flatten(), R.flatten()),
        values=Z.flatten()
    )
    
    max_density_points = []
    
    def calculate_marginal_pdf(value, dist_info):
        """Internal helper."""
        dist_name = dist_info['distribution']
        params = dist_info['params']
        try:
            if dist_name == 'genpareto':
                threshold = dist_info.get('threshold', 0)
                if value <= threshold:
                    return 0.0
                return stats.genpareto.pdf(value - threshold, *params)
            return getattr(stats, dist_name).pdf(value, *params)
        except Exception:
            return 1.0
    
    for T in return_periods:
        # logging omitted
        
        fig, ax = plt.subplots(figsize=(8, 6))
        
        Z_numeric = np.array(Z, dtype=float)
        
        try:
            cs = ax.contour(S, R, Z_numeric, levels=[T])
            
            all_vertices = []
            
            try:
                if hasattr(cs, 'collections') and len(cs.collections) > 0:
                    contour_paths = cs.collections[0].get_paths()
                    for path in contour_paths:
                        vertices = path.vertices
                        if len(vertices) > 0:
                            sea_levels = vertices[:, 0]
                            rainfalls = vertices[:, 1]
                            points = np.column_stack((sea_levels, rainfalls))
                            all_vertices.append(points)
                else:
                    raise AttributeError("collections not available")
            except (AttributeError, IndexError):
                try:
                    if hasattr(cs, 'allsegs') and len(cs.allsegs) > 0:
                        level_segments = cs.allsegs[0]
                        for segment in level_segments:
                            if len(segment) > 0:
                                vertices = np.array(segment)
                                sea_levels = vertices[:, 0]
                                rainfalls = vertices[:, 1]
                                points = np.column_stack((sea_levels, rainfalls))
                                all_vertices.append(points)
                    else:
                        raise AttributeError("allsegs not available")
                except (AttributeError, IndexError):
                    try:
                        for collection in cs.collections:
                            for path in collection.get_paths():
                                vertices = path.vertices
                                if len(vertices) > 0:
                                    sea_levels = vertices[:, 0]
                                    rainfalls = vertices[:, 1]
                                    points = np.column_stack((sea_levels, rainfalls))
                                    all_vertices.append(points)
                    except Exception as e3:
                        # logging omitted
                        plt.close(fig)
                        continue
            
            plt.close(fig)
            
            if not all_vertices:
                # logging omitted
                continue
                
            all_points = np.vstack(all_vertices)
            
            pdf_values = np.zeros(len(all_points))
            
            for i, (s_val, r_val) in enumerate(all_points):
                rainfall_cdf = calculate_probability(r_val, fits[rainfall_col])
                sea_level_cdf = calculate_probability(s_val, fits[sea_level_col])
                
                u = np.clip(rainfall_cdf, 0.001, 0.999)
                v = np.clip(sea_level_cdf, 0.001, 0.999)
                
                if copula_model is None:
                    rainfall_pdf = calculate_marginal_pdf(r_val, fits[rainfall_col])
                    sea_level_pdf = calculate_marginal_pdf(s_val, fits[sea_level_col])
                    pdf_values[i] = rainfall_pdf * sea_level_pdf
                else:
                    try:
                        is_survival_clayton = (best_copula == 'survival_clayton' and 
                                               copula_results[best_copula].get('method') == 'manual_transform')
                        
                        if is_survival_clayton:
                            u_transformed = 1 - u
                            v_transformed = 1 - v
                            copula_pdf = copula_model.pdf([u_transformed, v_transformed])
                        else:
                            copula_pdf = copula_model.pdf([u, v])
                            
                        rainfall_pdf = calculate_marginal_pdf(r_val, fits[rainfall_col])
                        sea_level_pdf = calculate_marginal_pdf(s_val, fits[sea_level_col])
                        
                        pdf_values[i] = copula_pdf * rainfall_pdf * sea_level_pdf
                    except Exception as e:
                        # logging omitted
                        pdf_values[i] = 0
            
            if np.all(pdf_values == 0):
                # logging omitted
                continue
                
            max_idx = np.argmax(pdf_values)
            max_point = all_points[max_idx]
            max_density_sea_level, max_density_rainfall = max_point
            
            point_return_period = contour_func(max_density_sea_level, max_density_rainfall)
            
            max_density_points.append({
                'ReturnPeriod': T,
                'SeaLevel_m': max_density_sea_level,
                'Rainfall_mm': max_density_rainfall,
                'PDF_Value': pdf_values[max_idx],
                'Verified_ReturnPeriod': float(point_return_period) if point_return_period is not None else np.nan
            })
            
            # logging omitted
            
        except Exception as e:
            # logging omitted
            plt.close(fig)
            continue
    
    if max_density_points:
        results_df = pd.DataFrame(max_density_points)
        # logging omitted
        return results_df
    else:
        # logging omitted
        return pd.DataFrame()

def plot_return_periods_single_station(ax, R, S, Z, compound_data, max_density_df,
                                      rainfall_threshold, sea_level_threshold, 
                                      station_name, title, x_max, y_max):
    """Internal helper."""
    # Define uniform color for points and background
    uniform_color = (80/255, 29/255, 138/255)  # Deep purple
    contour_color = '#46717C'  # Blue-green color for contours
    
    # Add background color matching the points but with high transparency
    ax.add_patch(plt.Rectangle((sea_level_threshold, rainfall_threshold), 
                              x_max - sea_level_threshold, y_max - rainfall_threshold, 
                              facecolor=uniform_color, alpha=0.15, zorder=0))
    
    # Set contour levels
    levels = [5, 10, 20, 50, 100, 200, 500, 1000]
    
    # Draw contours with increased thickness
    cs = ax.contour(S, R, Z, levels=levels, colors=contour_color, linewidths=2.0, 
                   corner_mask=False, antialiased=True)
    
    # Add contour labels with larger font
    labels = ax.clabel(cs, inline=False, fontsize=12, fmt='%d', colors=contour_color)
    
    # Adjust label positions
    for label in labels:
        pos = label.get_position()
        label.set_position((pos[0], pos[1] + 8))
    
    # Draw background points and extreme event points
    background_points = compound_data[~compound_data['Is_Compound_Extreme']]
    ax.scatter(
        background_points['Sea_Level'],
        background_points['Rainfall'],
        c='lightgray',
        s=15,
        alpha=0.5,
        edgecolors='gray',
        linewidths=0.5,
    )
    
    # Select compound extreme events
    compound_extreme = compound_data[compound_data['Is_Rep_Compound_Extreme'] == True]
    
    # Plot compound extreme events with larger size
    ax.scatter(
        compound_extreme['Sea_Level'],
        compound_extreme['Rainfall'],
        c=[uniform_color],
        s=45,
        alpha=0.7,
        edgecolors=None,
        linewidths=1.0,
        zorder=10,
    )
    
    if max_density_df is not None and not max_density_df.empty:
        triangle_color = '#e7298a'
        ax.scatter(
            max_density_df['SeaLevel_m'],
            max_density_df['Rainfall_mm'],
            c=triangle_color,
            marker='^',
            s=70,
            edgecolors=triangle_color,
            linewidths=1.5,
            zorder=30,
        )
    
    ax.set_xlim(sea_level_threshold, x_max)
    ax.set_ylim(rainfall_threshold, y_max)
    
    # Set x-axis ticks to one decimal place with increased font size
    import matplotlib.ticker as ticker
    ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%.1f'))
    
    # Set axis labels and title with larger font sizes
    ax.set_xlabel('Maximum Water Level (m)', fontsize=14, fontweight='bold')
    ax.set_ylabel('24 Hours Cumulative Rainfall (mm)', fontsize=14, fontweight='bold')
    ax.set_title(title, fontsize=16, fontweight='bold')
    
    # Add grid with increased alpha
    ax.grid(True, linestyle='--', alpha=0.4, color='gray', linewidth=0.8)
    
    # Make the borders (spines) much thicker
    for spine in ax.spines.values():
        spine.set_linewidth(2.0)
    
    # Adjust tick parameters - make them larger and thicker
    ax.tick_params(axis='both', which='major', width=2.0, length=6, labelsize=12, 
                  direction='out', pad=6)
    ax.tick_params(axis='both', which='minor', width=1.5, length=4, 
                  direction='out')
    
    # Set tick label font weight
    for tick in ax.get_xticklabels():
        tick.set_fontweight('bold')
    for tick in ax.get_yticklabels():
        tick.set_fontweight('bold')

def fit_marginal_distributions_selective(compound_data, column, min_samples=20, extreme_threshold=None,
                                        excluded_distributions=None):
    """Internal helper."""
    extreme_data_series = compound_data[compound_data['Is_Rep_Compound_Extreme'] == True][column].dropna()
    
    if len(extreme_data_series) < min_samples:
        if column == 'Rainfall':
            extra_data = compound_data[(compound_data['Is_Extreme_Rainfall'] == True) & 
                                     (compound_data['Is_Rep_Compound_Extreme'] == False)][column].dropna()
            extreme_data_series = pd.concat([extreme_data_series, extra_data.sample(min(len(extra_data), min_samples - len(extreme_data_series)), random_state=1)])
        else:  # Sea_Level
            extra_data = compound_data[(compound_data['Is_Extreme_Sea_Level'] == True) & 
                                     (compound_data['Is_Rep_Compound_Extreme'] == False)][column].dropna()
            extreme_data_series = pd.concat([extreme_data_series, extra_data.sample(min(len(extra_data), min_samples - len(extreme_data_series)), random_state=1)])
        extreme_data_series = extreme_data_series.dropna()

    if len(extreme_data_series) < 10:
        loc, scale = np.mean(extreme_data_series) if len(extreme_data_series) > 0 else 0, np.std(extreme_data_series) if len(extreme_data_series) > 1 else 1.0
        return {
            'distribution': 'norm', 'params': (loc, scale), 'bic': float('inf'),
            'ks_stat': float('inf'), 'ks_pvalue': 0.0,
            'num_samples_for_pot_calc_prob': len(extreme_data_series),
            'exceedance_values_for_calc_prob': np.array([])
        }

    extreme_data_np_array = extreme_data_series.values
    n_samples = len(extreme_data_np_array)

    gpd_actual_threshold = extreme_threshold
    if gpd_actual_threshold is None:
        if column == 'Rainfall':
            gpd_actual_threshold = np.percentile(extreme_data_np_array, 10) if len(extreme_data_np_array) > 0 else 0
        else: 
            gpd_actual_threshold = np.min(extreme_data_np_array) if len(extreme_data_np_array) > 0 else 0
    
    fit_results = []
    
    def calculate_bic(k, loglik, n):
        return k * np.log(n) - 2 * loglik

    common_props_for_prob_calc = {
        'num_samples_for_pot_calc_prob': n_samples
    }

    candidate_distributions = resolve_marginal_candidates(excluded_distributions)

    for dist_name in candidate_distributions:
        try:
            if dist_name == 'genpareto':
                exceedances = extreme_data_np_array[extreme_data_np_array > gpd_actual_threshold] - gpd_actual_threshold
                exceedances_for_fit = exceedances[exceedances > 1e-9]
                if len(exceedances_for_fit) < 5:
                    continue
                params = stats.genpareto.fit(exceedances_for_fit, floc=0)
                loglik = np.sum(stats.genpareto.logpdf(exceedances_for_fit, *params))
                bic = calculate_bic(len(params)-1 if len(params)>1 else 1, loglik, len(exceedances_for_fit))
                ks_stat, ks_pvalue = stats.kstest(exceedances_for_fit, lambda x: stats.genpareto.cdf(x, *params))
                cm_res = stats.cramervonmises(exceedances_for_fit, 'genpareto', args=params)
                fit_results.append({
                    'distribution': 'genpareto', 'params': params, 'loglik': loglik, 'bic': bic,
                    'ks_stat': ks_stat, 'ks_pvalue': ks_pvalue,
                    'cm_stat': cm_res.statistic, 'cm_pvalue': cm_res.pvalue,
                    'threshold': gpd_actual_threshold,
                    'exceedance_values_for_calc_prob': exceedances_for_fit,
                    **common_props_for_prob_calc
                })
            else:
                dist_obj = getattr(stats, dist_name)
                params = dist_obj.fit(extreme_data_np_array)
                loglik = np.sum(dist_obj.logpdf(extreme_data_np_array, *params))
                bic = calculate_bic(len(params), loglik, n_samples)
                ks_stat, ks_pvalue = stats.kstest(extreme_data_np_array, lambda x: dist_obj.cdf(x, *params))
                cm_res = stats.cramervonmises(extreme_data_np_array, dist_name, args=params)
                fit_results.append({**{'distribution': dist_name, 'params': params, 'loglik': loglik, 'bic': bic,
                    'ks_stat': ks_stat, 'ks_pvalue': ks_pvalue,
                    'cm_stat': cm_res.statistic, 'cm_pvalue': cm_res.pvalue},
                    **common_props_for_prob_calc})
        except Exception:
            continue

    if not fit_results:
        # logging omitted
        loc, scale = np.mean(extreme_data_np_array) if len(extreme_data_np_array) > 0 else 0, np.std(extreme_data_np_array) if len(extreme_data_np_array) > 1 else 1.0
        return {
            'distribution': 'norm', 'params': (loc, scale), 'bic': float('inf'),
            'ks_stat': float('inf'), 'ks_pvalue': 0.0, 
            **common_props_for_prob_calc,
            'exceedance_values_for_calc_prob': np.array([])
        }

    ks_significance_level = 0.05
    valid_fits_by_ks = [fit for fit in fit_results if fit.get('ks_pvalue', 0) > ks_significance_level]
    selected_fits_for_bic = valid_fits_by_ks if valid_fits_by_ks else fit_results
    best_fit = min(selected_fits_for_bic, key=lambda x: x['bic'])
    
    # logging omitted
    return best_fit

def fit_copula_models_selective(compound_data, rainfall_col, sea_level_col, rainfall_threshold=None, sea_level_threshold=None,
                               excluded_distributions=None):
    """Internal helper."""
    # Calculate 95th percentile if thresholds not provided
    if rainfall_threshold is None:
        rainfall_threshold = np.percentile(compound_data[rainfall_col], 95)
    if sea_level_threshold is None:
        sea_level_threshold = np.percentile(compound_data[sea_level_col], 95)
        
    print(f"Rainfall threshold: {rainfall_threshold:.2f} mm")
    print(f"Sea level threshold: {sea_level_threshold:.2f} m")
    
    # Select compound extreme events
    extreme_events = compound_data[compound_data['Is_Rep_Compound_Extreme'] == True].copy()
    
    print(f"Compound extreme events: {len(extreme_events)}")
    
    # Check if enough extreme events
    if len(extreme_events) < 20:
        print(f"Warning: Few compound extreme events ({len(extreme_events)}), will use enhanced fitting strategy")
    
    print(f"Using {len(extreme_events)} compound extreme events to fit Copula models")
    
    fits = {
        rainfall_col: fit_marginal_distributions_selective(compound_data, rainfall_col, extreme_threshold=rainfall_threshold, excluded_distributions=excluded_distributions),
        sea_level_col: fit_marginal_distributions_selective(compound_data, sea_level_col, extreme_threshold=sea_level_threshold, excluded_distributions=excluded_distributions)
    }
    
    print("\nMarginal distribution fitting results (Best by K-S & BIC):")
    for col, fit in fits.items():
        dist_name = fit['distribution']
        params = fit['params']
        bic_val = fit['bic']
        
        param_str_parts = []
        if dist_name == 'genpareto':
            shape, _, scale = params
            threshold = fit['threshold']
            param_str_parts.extend([f"shape={shape:.4f}", f"scale={scale:.4f}", f"threshold={threshold:.4f}"])
        elif dist_name == 'norm':
            loc, scale = params
            param_str_parts.extend([f"loc={loc:.4f}", f"scale={scale:.4f}"])
        elif dist_name == 'genextreme':
            shape, loc, scale = params
            param_str_parts.extend([f"shape={shape:.4f}", f"loc={loc:.4f}", f"scale={scale:.4f}"])
        elif dist_name == 'gamma':
            shape, loc, scale = params
            param_str_parts.extend([f"shape={shape:.4f}", f"loc={loc:.4f}", f"scale={scale:.4f}"])
        elif dist_name == 'lognorm':
            shape, loc, scale = params
            param_str_parts.extend([f"shape={shape:.4f}", f"loc={loc:.4f}", f"scale={scale:.4f}"])
        elif dist_name == 'logistic':
            loc, scale = params
            param_str_parts.extend([f"loc={loc:.4f}", f"scale={scale:.4f}"])
        elif dist_name == 'gumbel_r':
            loc, scale = params
            param_str_parts.extend([f"loc={loc:.4f}", f"scale={scale:.4f}"])
        else:
            param_str_parts.append(f"params={params}")
        
        ks_p_val_str = f"{fit.get('ks_pvalue', -1):.4f}" if 'ks_pvalue' in fit else "N/A"
        cm_p_val_str = f"{fit.get('cm_pvalue', -1):.4f}" if 'cm_pvalue' in fit else "N/A"
        print(f"  {col}: {dist_name.capitalize()}, BIC={bic_val:.2f}, K-S p-val={ks_p_val_str}, CM p-val={cm_p_val_str}, Params: {', '.join(param_str_parts)}")

    # Transform data using theoretical CDFs from fitted marginal distributions
    print("\nTransforming data to [0,1] using theoretical CDFs...")
    u_values = extreme_events[rainfall_col].apply(lambda x: calculate_probability(x, fits[rainfall_col]))
    v_values = extreme_events[sea_level_col].apply(lambda x: calculate_probability(x, fits[sea_level_col]))
    
    # Prepare data - ensure strictly in (0,1) interval
    u = np.clip(u_values.values, 0.001, 0.999)
    v = np.clip(v_values.values, 0.001, 0.999)
    data = np.column_stack((u, v))

    # Diagnostic output for transformed data
    print(f"Transformed Rainfall (u) range: [{u.min():.6f}, {u.max():.6f}]")
    print(f"Transformed Sea Level (v) range: [{v.min():.6f}, {v.max():.6f}]")

    results = {}
    
    # Helper function to calculate NLL using PDF
    def calculate_nll_from_pdf(copula_model, data, model_name):
        log_lik = 0
        pdf_values = np.zeros(len(data))
        if not hasattr(copula_model, 'pdf'):
             raise AttributeError(f"Critical Error: pdf method not found for {model_name} Copula. Program will terminate.")

        for i in range(len(data)):
             try:
                  pdf_val = copula_model.pdf(data[i])
                  if not np.isfinite(pdf_val) or pdf_val <= 0:
                      raise ValueError(f"Critical Error: PDF for {model_name} at point {i} is invalid (value: {pdf_val}). Program will terminate.")
                  pdf_values[i] = pdf_val
             except Exception as pdf_err:
                  raise RuntimeError(f"Critical Error during PDF calculation for {model_name} at point {i}: {pdf_err}. Program will terminate.") from pdf_err
        
        try:
            if np.any(pdf_values <= 0):
                 raise ValueError(f"Critical Error: PDF values for {model_name} include non-positive numbers before log. Program will terminate.")
            log_pdf_values = np.log(pdf_values)
            if np.any(np.isinf(log_pdf_values)) or np.any(np.isnan(log_pdf_values)):
                raise ValueError(f"Critical Error: Log of PDF for {model_name} resulted in Inf or NaN. Program will terminate.")
            log_lik = np.sum(log_pdf_values)
        except ValueError as ve:
            raise ve
        except Exception as log_err:
            raise RuntimeError(f"Critical Error during log-likelihood sum for {model_name}: {log_err}. Program will terminate.") from log_err
        
        nll = -log_lik
        if not np.isfinite(nll):
            raise ValueError(f"Critical Error: NLL for {model_name} is not finite (value: {nll}). Program will terminate.")
        return nll

    try:
        # Gaussian Copula
        cop_gaussian = GaussianCopula(dim=2)
        cop_gaussian.fit(data)
        nll = calculate_nll_from_pdf(cop_gaussian, data, "Gaussian")
        k = len(np.ravel(cop_gaussian.params))
        n = len(data)
        bic = k * np.log(n) + 2 * nll
        cm_p_value, cm_stat = perform_cramer_von_mises_test(data, cop_gaussian)
        ks_p_value, ks_stat = perform_kolmogorov_smirnov_test(data, cop_gaussian)
        results['gaussian'] = {
            'params': cop_gaussian.params, 'bic': bic, 'nll': nll, 'model': cop_gaussian,
            'cm_p_value': cm_p_value, 'cm_stat': cm_stat, 'ks_p_value': ks_p_value, 'ks_stat': ks_stat
        }
        
        # Student-t Copula
        cop_student = StudentCopula(dim=2)
        cop_student.fit(data)
        if isinstance(cop_student.params, tuple):
            rho_matrix, df = cop_student.params
            rho = rho_matrix[0, 1] if rho_matrix.ndim == 2 else rho_matrix
        else:
            rho = cop_student.params[0, 1]
            df = getattr(cop_student, 'df', 3.0)
        rho = np.clip(float(rho), -0.999, 0.999)
        df = max(float(df), 2.01)
        nll = calculate_nll_from_pdf(cop_student, data, "Student-t")
        k = 2
        bic = k * np.log(n) + 2 * nll
        cm_p_value, cm_stat = perform_cramer_von_mises_test(data, cop_student)
        ks_p_value, ks_stat = perform_kolmogorov_smirnov_test(data, cop_student)
        results['t'] = {
            'params': (rho, df), 'bic': bic, 'nll': nll, 'model': cop_student,
            'cm_p_value': cm_p_value, 'cm_stat': cm_stat, 'ks_p_value': ks_p_value, 'ks_stat': ks_stat
        }
        
        # Clayton Copula
        cop_clayton = ClaytonCopula(dim=2)
        cop_clayton.fit(data)
        nll = calculate_nll_from_pdf(cop_clayton, data, "Clayton")
        k = 1
        bic = k * np.log(n) + 2 * nll
        cm_p_value, cm_stat = perform_cramer_von_mises_test(data, cop_clayton)
        ks_p_value, ks_stat = perform_kolmogorov_smirnov_test(data, cop_clayton)
        results['clayton'] = {
            'params': cop_clayton.params, 'bic': bic, 'nll': nll, 'model': cop_clayton,
            'cm_p_value': cm_p_value, 'cm_stat': cm_stat, 'ks_p_value': ks_p_value, 'ks_stat': ks_stat
        }
        
        # Gumbel Copula
        cop_gumbel = GumbelCopula(dim=2)
        cop_gumbel.fit(data)
        nll = calculate_nll_from_pdf(cop_gumbel, data, "Gumbel")
        k = 1
        bic = k * np.log(n) + 2 * nll
        cm_p_value, cm_stat = perform_cramer_von_mises_test(data, cop_gumbel)
        ks_p_value, ks_stat = perform_kolmogorov_smirnov_test(data, cop_gumbel)
        results['gumbel'] = {
            'params': cop_gumbel.params, 'bic': bic, 'nll': nll, 'model': cop_gumbel,
            'cm_p_value': cm_p_value, 'cm_stat': cm_stat, 'ks_p_value': ks_p_value, 'ks_stat': ks_stat
        }
        
        # Frank Copula
        cop_frank = FrankCopula(dim=2)
        cop_frank.fit(data)
        nll = calculate_nll_from_pdf(cop_frank, data, "Frank")
        k = 1
        bic = k * np.log(n) + 2 * nll
        cm_p_value, cm_stat = perform_cramer_von_mises_test(data, cop_frank)
        ks_p_value, ks_stat = perform_kolmogorov_smirnov_test(data, cop_frank)
        results['frank'] = {
            'params': cop_frank.params, 'bic': bic, 'nll': nll, 'model': cop_frank,
            'cm_p_value': cm_p_value, 'cm_stat': cm_stat, 'ks_p_value': ks_p_value, 'ks_stat': ks_stat
        }
        
        # Joe Copula
        cop_joe = SimpleJoeCopula(dim=2)
        cop_joe.fit(data)
        nll = calculate_nll_from_pdf(cop_joe, data, "Joe")
        k = 1
        bic = k * np.log(n) + 2 * nll
        cm_p_value, cm_stat = perform_cramer_von_mises_test(data, cop_joe)
        ks_p_value, ks_stat = perform_kolmogorov_smirnov_test(data, cop_joe)
        results['joe'] = {
            'params': cop_joe.params, 'bic': bic, 'nll': nll, 'model': cop_joe,
            'cm_p_value': cm_p_value, 'cm_stat': cm_stat, 'ks_p_value': ks_p_value, 'ks_stat': ks_stat
        }
        
        # Survival Clayton Copula
        u_orig = data[:, 0]
        v_orig = data[:, 1]
        u_survival_transform = 1 - u_orig
        v_survival_transform = 1 - v_orig
        data_for_survival_fit = np.column_stack((u_survival_transform, v_survival_transform))
        cop_clayton_std_for_survival = ClaytonCopula(dim=2)
        cop_clayton_std_for_survival.fit(data_for_survival_fit)
        theta_sc = cop_clayton_std_for_survival.params
        nll = calculate_nll_from_pdf(cop_clayton_std_for_survival, data_for_survival_fit, "Survival Clayton (manual)")
        k = 1
        bic = k * np.log(n) + 2 * nll
        cm_p_value, cm_stat = perform_cramer_von_mises_test(data_for_survival_fit, cop_clayton_std_for_survival)
        ks_p_value, ks_stat = perform_kolmogorov_smirnov_test(data_for_survival_fit, cop_clayton_std_for_survival)
        results['survival_clayton'] = {
            'params': theta_sc, 'bic': bic, 'nll': nll, 'model': cop_clayton_std_for_survival, 
            'cm_p_value': cm_p_value, 'cm_stat': cm_stat, 'ks_p_value': ks_p_value, 'ks_stat': ks_stat,
            'method': 'manual_transform'
        }
        
    except Exception as e:
        print(f"Error in Copula fitting: {e}")
        results['independence'] = {
            'params': [0], 'bic': float('inf'), 'nll': float('inf'), 'model': None,
            'cm_p_value': 1.0, 'cm_stat': 0.0, 'ks_p_value': 1.0, 'ks_stat': 0.0
        }
        return results, 'independence', fits

    # Model selection logic
    print("\n--- Copula Model Evaluation Summary ---")
    for name, res_dict in results.items():
        param_str = f"{res_dict['params']}"
        if name == 't':
             param_str = f"rho={res_dict['params'][0]:.3f}, df={res_dict['params'][1]:.1f}"
        elif name == 'gaussian':
             if isinstance(res_dict['params'], np.ndarray) and res_dict['params'].ndim == 2:
                 param_str = f"corr={res_dict['params'][0,1]:.4f}" 
             else:
                 param_str = f"{res_dict['params']}" 
        elif name == 'survival_clayton':
            param_str = f"theta={res_dict['params']:.4f}"
        else:
             param_str = f"{res_dict['params']:.4f}"
        print(f"  {name.capitalize()} Copula: BIC={res_dict['bic']:.2f}, K-S p-val={res_dict.get('ks_p_value', float('nan')):.4f}, Params: {param_str}")

    # Model selection
    ks_significance_level = 0.05
    analyzable_models = {
        name: res for name, res in results.items()
        if pd.notna(res.get('bic')) and pd.notna(res.get('ks_p_value'))
    }

    if not analyzable_models:
        raise RuntimeError("CRITICAL: No Copula models have valid BIC and K-S p-value for selection. Program will terminate.")

    passed_ks_test_models = {
        name: res for name, res in analyzable_models.items()
        if res['ks_p_value'] > ks_significance_level
    }

    if passed_ks_test_models:
        print(f"\nModels passing K-S test (p > {ks_significance_level}): {list(passed_ks_test_models.keys())}")
        best_copula_name = min(passed_ks_test_models.items(), key=lambda x: x[1]['bic'])[0]
        print(f"Selected best model (passes K-S test & lowest BIC): {best_copula_name.capitalize()}")
    else:
        print(f"\nWarning: No Copula model passed the K-S test (p > {ks_significance_level}).")
        best_copula_name = min(analyzable_models.items(), key=lambda x: x[1]['bic'])[0]
        print(f"Selected best model based on lowest BIC (as fallback, K-S test not passed): {best_copula_name.capitalize()}")

    return results, best_copula_name, fits

def test_datum_correction(tide_file_path, offset_file_path):
    """Internal helper."""
    # logging omitted
    
    # logging omitted
    try:
        original_data = pd.read_csv(tide_file_path, encoding='utf-8')
        # logging omitted
        
        date_cols_tide = ['Year', 'Month', 'Day']
        for col in date_cols_tide:
            original_data[col] = pd.to_numeric(original_data[col], errors='coerce')
        original_data = original_data.dropna(subset=date_cols_tide)
        for col in date_cols_tide:
            original_data[col] = original_data[col].astype(int)
        
        sea_level_col_original = 'Sea level(m)'
        sea_level_col_target = 'Sea_Level'
        if sea_level_col_original in original_data.columns:
            original_data.rename(columns={sea_level_col_original: sea_level_col_target}, inplace=True)
        
        original_data[sea_level_col_target] = pd.to_numeric(original_data[sea_level_col_target], errors='coerce')
        original_data = original_data.dropna(subset=[sea_level_col_target])
        
        
        # logging omitted
        # logging omitted
        # logging omitted
        # logging omitted
        # logging omitted
        # logging omitted
        
    except Exception as e:
        # logging omitted
        return None, None
    
    # logging omitted
    offset_data = load_datum_correction_file(offset_file_path)
    if offset_data is None:
        # logging omitted
        return original_data, None
    
    # logging omitted
    corrected_data = apply_datum_correction(original_data, offset_data, sea_level_col_target)
    
    # logging omitted
    if f'{sea_level_col_target}_Original' in corrected_data.columns:
        original_stats = corrected_data[f'{sea_level_col_target}_Original'].describe()
        corrected_stats = corrected_data[sea_level_col_target].describe()
        
        # logging omitted
        # logging omitted
        # logging omitted
        # logging omitted
        # logging omitted
        
        # logging omitted
        # logging omitted
        # logging omitted
        # logging omitted
        # logging omitted
        
        mean_correction = corrected_stats['mean'] - original_stats['mean']
        # logging omitted
        # logging omitted
        # logging omitted
    
    if corrected_data is not None:
        output_file = tide_file_path.replace('.csv', '_datum_corrected.csv')
        corrected_data.to_csv(output_file, index=False, encoding='utf-8-sig')
        # logging omitted
    
    # logging omitted
    return original_data, corrected_data

def main_multi_station():    
    """Internal helper."""
    station_configs = [
        {
            'name': 'HD06', 
            'tide_file': str(DATA_DIR / 'hd06_lat35.56805_lon140.04555_tide.csv'),
            'rainfall_files': [str(DATA_DIR / '196401_202412_rainfall_47682.csv')],
            'offset_file': str(DATA_DIR / 'hd006_autum.csv'),  # Vertical datum offset file
            'title': 'HD06',
            'x_max': 4.0,  # Plot x-axis limit (m)
            'y_max': 350,  # Plot y-axis limit (mm)
            'excluded_distributions': ['gamma', 'genpareto'],
        },
        {
            'name': 'MA15',
            'tide_file': str(DATA_DIR / 'ma15_lat35.64862_lon139.77_tide.csv'),
            'rainfall_files': [str(DATA_DIR / '196401_202412_rainfall_47662.csv'),str(DATA_DIR / '197601_202412_rainfall_0370.csv')],
            'offset_file': str(DATA_DIR / 'ma15_autum.csv'),  # Vertical datum offset file
            'title': 'MA15',
            'x_max': 2.5,  # Plot x-axis limit (m)
            'y_max': 300,  # Plot y-axis limit (mm)
        },
        {
            'name': 'HD08',
            'tide_file': str(DATA_DIR / 'hd08_lat35.28805_lon139.65138_tide.csv'),
            'rainfall_files': [str(DATA_DIR / '197601_202412_rainfall_0392.csv'),str(DATA_DIR / '197601_202412_rainfall_0391.csv')],
            'offset_file': str(DATA_DIR / 'hd08_autum.csv'),  # Vertical datum offset file
            'title': 'HD08',
            'x_max': 3.0,  # Plot x-axis limit (m)
            'y_max': 250,  # Plot y-axis limit (mm)
            'excluded_distributions': ['gamma', 'genpareto'],
        },
        {
            'name': 'GS01',
            'tide_file': str(DATA_DIR / 'gs01_lat35.16_lon139.6155_tide.csv'),
            'rainfall_files': [str(DATA_DIR / '197601_202412_rainfall_0392.csv')],
            'offset_file': str(DATA_DIR / 'gs01_autum.csv'),  # Vertical datum offset file
            'title': 'GS01',
            'x_max': 2.0,  # Plot x-axis limit (m)
            'y_max': 300,  # Plot y-axis limit (mm)
        },
        {
            'name': 'MA14',
            'tide_file': str(DATA_DIR / 'ma14_lat34.91888_lon139.825_tide.csv'),
            'rainfall_files': [str(DATA_DIR / '196401_202412_rainfall_47672.csv')],
            'offset_file': str(DATA_DIR / 'ma14_autum.csv'),  # Vertical datum offset file
            'title': 'MA14',
            'x_max': 3.0,  # Plot x-axis limit (m)
            'y_max': 400,  # Plot y-axis limit (mm)
        },
        {
            'name': 'MA59',
            'tide_file': str(DATA_DIR / 'ma59_lat35.5_lon139.76667_tide.csv'),
            'rainfall_files': [str(DATA_DIR / '196401_202412_rainfall_0371.csv'),str(DATA_DIR / '196401_202412_rainfall_47670.csv')],
            'offset_file': str(DATA_DIR / 'ma59_autum.csv'),  # Vertical datum offset file
            'title': 'MA59',
            'x_max': 2.0,  # Plot x-axis limit (m)
            'y_max': 350,  # Plot y-axis limit (mm)
        }
    ]
    
    # logging omitted
    
    fig = create_multi_station_figure(station_configs, figsize=(18, 14))
    
    # logging omitted
    return fig

def demo_datum_correction():
    """Internal helper."""
    # logging omitted
    
    tide_file = str(DATA_DIR / 'ma15_lat35.64862_lon139.77_tide.csv')
    offset_file = str(DATA_DIR / 'ma15_autum.csv')
    
    # logging omitted
    # logging omitted
    # logging omitted
    
    original_data, corrected_data = test_datum_correction(tide_file, offset_file)
    
    if corrected_data is not None and 'Applied_Offset' in corrected_data.columns:
        applied_corrections = (corrected_data['Applied_Offset'] != 0).sum()
        if applied_corrections > 0:
            export_corrected_tide_data(corrected_data, 'MA15')
            # logging omitted
    
    if original_data is not None and corrected_data is not None:
        # logging omitted
        # logging omitted
        pass
    else:
        # logging omitted
        # logging omitted
    
        pass
    return original_data, corrected_data

def demo_contour_intersections():
    """Internal helper."""
    # logging omitted
    
    station_name = 'HD06'
    tide_file = str(DATA_DIR / 'hd06_lat35.56805_lon140.04555_tide.csv')
    rainfall_files = [str(DATA_DIR / '196401_202412_rainfall_47682.csv')]
    offset_file = str(DATA_DIR / 'hd006_autum.csv')
    
    # logging omitted
    # logging omitted
    # logging omitted
    # logging omitted
    
    try:
        tide_data, rainfall_data = load_data(tide_file, rainfall_files, export_rainfall=False, 
                                           offset_file_path=offset_file)
        
        if tide_data is None or rainfall_data is None:
            # logging omitted
            return None
        
        compound_data, rainfall_threshold, sea_level_threshold = create_compound_events(
            tide_data, rainfall_data, 
            rainy_day_def=0.1,
            rainfall_threshold_method='all_days'
        )
        
        if compound_data is None:
            # logging omitted
            return None
        
        copula_results, best_copula, fits = fit_copula_models_selective(
            compound_data, 'Rainfall', 'Sea_Level', 
            rainfall_threshold=rainfall_threshold, 
            sea_level_threshold=sea_level_threshold,
            excluded_distributions=['gamma', 'genpareto']
        )
        
        x_max, y_max = 4.0, 350
        R, S, Z = calculate_contours(
            compound_data, 'Rainfall', 'Sea_Level',
            copula_results, best_copula, fits,
            rainfall_threshold=rainfall_threshold,
            sea_level_threshold=sea_level_threshold,
            x_min=sea_level_threshold, x_max=x_max,
            y_min=rainfall_threshold, y_max=y_max,
            npoints=200
        )
        
        intersections_df = calculate_contour_axis_intersections(
            R, S, Z, 
            return_periods=[5, 10, 20, 50, 100, 200, 500, 1000]
        )
        
        if intersections_df is not None and not intersections_df.empty:
            intersections_df['Station'] = station_name
            intersections_df = intersections_df[['Station', 'ReturnPeriod_years', 'SeaLevel_at_zero_rainfall_m', 'Rainfall_at_zero_sealevel_mm', 'Notes']]
        
        if intersections_df is not None and not intersections_df.empty:
            demo_filename = f'demo_contour_intersections_{station_name}.csv'
            intersections_df.to_csv(demo_filename, index=False, encoding='utf-8-sig')
            # logging omitted
            # logging omitted
            # logging omitted
            # logging omitted
            print(intersections_df.head(10))
            # logging omitted
            axis_counts = intersections_df['Axis'].value_counts()
            for axis, count in axis_counts.items():
                # logging omitted
                pass
            return intersections_df
        else:
            # logging omitted
            # logging omitted
            return None
            
    except Exception as e:
        # logging omitted
        # logging omitted
        import traceback
        traceback.print_exc()
        return None

if __name__ == "__main__":
    import sys
    
    if len(sys.argv) > 1 and sys.argv[1] == '--test-datum':
        demo_datum_correction()
    elif len(sys.argv) > 1 and sys.argv[1] == '--test-intersections':
        demo_contour_intersections()
    else:
        multi_fig = main_multi_station()
        
        plt.show()