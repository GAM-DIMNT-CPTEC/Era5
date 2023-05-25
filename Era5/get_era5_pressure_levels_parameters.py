#! /usr/bin/env python3

import cdsapi

# Usar o ambiente CDSAPI

c = cdsapi.Client()

c.retrieve(
    'reanalysis-era5-pressure-levels',
    {
        'product_type': 'reanalysis',
        'format': 'grib',
        'variable': [
            'u_component_of_wind', 'v_component_of_wind', 'temperature',
            'specific_humidity', 'geopotential',
        ],
        'year': '2023',
        'month': [
            '02', '03',
        ],
        'day': [
            '01', '02', '03',
            '04', '05', '06',
            '07', '08', '09',
            '10', '11', '12',
            '13', '14', '15',
            '16', '17', '18',
            '19', '20', '21',
            '22', '23', '24',
            '25', '26', '27',
            '28', '29', '30',
            '31',
        ],
        'time': [
            '00:00', '12:00',
        ],
        'pressure_level': [
            '1000', 
            '925', 
            '850',  
            '500', 
            '250',
        ]
    },
    'Era520230201002023033112.atm.grib')
