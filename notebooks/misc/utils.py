# http://stackoverflow.com/questions/1094841/reusable-library-to-get-human-readable-version-of-file-size
def sizeof_fmt(num, sep=''):
    '''
    :param sep': separate between the number of unit, sometimes it's preferred
    for better readability
    '''
    for unit in ['Bytes', 'KB', 'MB', 'GB', 'TB', 'PB', 'EB', 'ZB']:
        if abs(num) < 1024.0:
            return "%3.1f %s" % (num, unit)
        num /= 1024.0
    return "%.1f %s" % (num, 'YB')


def count(df, colname):
    '''return the count and percentage for a given column after groupby'''
    # it can be any column other than 'study'
    if isinstance(colname, list):
        cols = colname + ['study']
    else:
        cols = ['study', colname]
    res = df[cols].groupby(colname).count().sort_values('study')
    res['percent'] = (res.study / res.sum().values[0]).apply('{0:.2%}'.format)
    return res
