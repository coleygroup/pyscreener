def time_stats(total_time, n):
    m, s = divmod(int(total_time), 60)
    h, m = divmod(m, 60)
    avg = total_time / n
    
    return h, m, s, avg