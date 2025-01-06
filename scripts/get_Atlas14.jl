"""web scrape Atlas 14 estimates"""

# Function to extract and parse the quantiles from the fetched data
function parse_quantiles_from_url(url)
    # Fetching data from the URL
    response = HTTP.get(url, require_ssl_verification=false)
    data = String(response.body)

    # Extracting the quantiles data
    start_idx = findfirst("quantiles = ", data)
    end_idx = findnext(";", start_idx) - 1
    quantiles_str = data[start_idx:end_idx]

    # Cleaning up the quantiles string for JSON parsing
    quantiles_str = replace(quantiles_str, "quantiles = " => "")
    quantiles_str = replace(quantiles_str, "'" => "\"")

    # Parsing the JSON string to a Julia array
    quantiles_array = JSON.parse(quantiles_str)

    # Converting the string values to floats
    quantiles_floats = map(x -> parse(Float64, x), quantiles_array)

    return quantiles_floats
end

function get_Atlas14(lat, lon, row, col)
    # row is linked with durations
    url = "https://hdsc.nws.noaa.gov/cgi-bin/hdsc/new/cgi_readH5.py?lat=" * lat * "&lon=" * lon * "&type=pf&data=depth&units=english&series=pds"
    a = DataFrame(CSV.File(download(url)))[1, 1]
    a = a[15:length(a)-2]
    a = replace(a, "'" => "")
    b = split(a, "], [")
    c1 = [parse(Float64, s) for s in split(b[row], ", ")][col]
    return c1
end

function get_Atlas14_IDF(lat, lon)
    url = "https://hdsc.nws.noaa.gov/cgi-bin/hdsc/new/cgi_readH5.py?lat=" * lat * "&lon=" * lon * "&type=pf&data=depth&units=english&series=pds"
    a = DataFrame(CSV.File(download(url)))[1, 1]
    a = a[15:length(a)-2]
    a = replace(a, "'" => "")
    b = split(a, "], [")
    rls = [[parse(Float64, s) for s in split(b[10], ", ")][i] for i in [1, 2, 3, 4, 5, 6, 7, 8, 9]]
    return rls
end

