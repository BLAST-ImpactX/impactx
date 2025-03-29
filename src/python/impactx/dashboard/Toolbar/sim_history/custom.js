window.getSimStatusColor = function(status) {
    if (status === 'Completed') return 'success';
    if (status === 'In Progress') return 'info';
    if (status === 'Failed') return 'error';
    return 'default';
}
