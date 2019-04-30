function update_db(filename)

    load(filename, 'mogs')
    for n=1:size(mogs,2)
        mog = Mog(mogs(n));
        mogs(n) = mog;
    end
    save(filename, 'mogs', '-append')
end
