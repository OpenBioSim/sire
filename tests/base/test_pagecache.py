import sire as sr


def test_pagecache(ala_mols, tmpdir):
    import pathlib

    mols = ala_mols

    PageCache = sr.base.PageCache

    root_dir = pathlib.PurePath(PageCache.root_directory()).as_posix()

    assert pathlib.Path(root_dir).exists()

    # should be the current path
    cwd = pathlib.Path().cwd().as_posix()

    assert root_dir == cwd

    d = tmpdir.mkdir("test_pagecache")

    new_root_dir = pathlib.PurePath(d.strpath).as_posix()

    PageCache.set_root_directory(new_root_dir)

    root_dir2 = pathlib.PurePath(PageCache.root_directory()).as_posix()

    assert pathlib.Path(root_dir2).exists()

    # should be the tempdir
    assert root_dir2 == new_root_dir

    # restore to the original root dir
    PageCache.set_root_directory(root_dir)

    assert pathlib.PurePath(PageCache.root_directory()).as_posix() == root_dir

    c = PageCache("temp_pagecache_XXXXXX", 32 * 1024)

    assert c.page_size() == 32 * 1024
    assert c.num_pages() == 0
    assert c.num_bytes() == 0

    assert "temp_pagecache_" in c.cache_dir()

    h1 = c.store(42)

    assert h1.fetch() == 42

    h2 = c.store("Hello Python Page Cache")

    assert h2.fetch() == "Hello Python Page Cache"

    assert h1.fetch() == 42

    h3 = c.store(mols)

    mols2 = h3.fetch()

    assert mols2.num_molecules() == mols.num_molecules()
    assert mols2.num_atoms() == mols.num_atoms()

    assert c.num_pages() > 0
    cache_dir = c.cache_dir()

    assert "temp_pagecache_" in cache_dir

    handles = []

    for i in range(0, 5000):
        handles.append(c.store(i))

    for i, h in enumerate(handles):
        assert h.fetch() == i

    # the cache may not have flushed to pages or disk yet
    # so we should wait a little bit of time to let it do that...
    import time

    n_sleeps = 0

    while c.num_pages() < 2:
        time.sleep(0.1)
        n_sleeps += 1

        if n_sleeps > 50:
            break

    assert c.num_pages() >= 2

    cache_dir = pathlib.PurePath(c.cache_dir()).as_posix()

    assert "temp_pagecache_" in cache_dir

    assert pathlib.Path(cache_dir).exists()
