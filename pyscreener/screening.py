from pyscreener import docking, md, dft

def screen(mode, **kwargs):
    if mode == 'docking':
        return docking.dock(**kwargs)
    if mode == 'md':
        return md.simulate(**kwargs)
    if mode == 'dft':
        return dft.calculate(**kwargs)
    
    raise ValueError