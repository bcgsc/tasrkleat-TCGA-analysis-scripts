# http://stackoverflow.com/questions/1094841/reusable-library-to-get-human-readable-version-of-file-size
def sizeof_fmt(num, sep=''):
    '''
    :param sep': separate between the number of unit, sometimes it's preferred
    for better readability
    '''
    for unit in ['Bytes','KB','MB','GB','TB','PB','EB','ZB']:
        if abs(num) < 1024.0:
            return "%3.1f %s" % (num, unit)
        num /= 1024.0
    return "%.1f %s" % (num, 'YB')
